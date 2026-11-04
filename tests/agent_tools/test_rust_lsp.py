"""Real subprocess/pipe protocol tests against a synthetic LSP, not rust-analyzer."""
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("rust_lsp", ROOT / "scripts/rust_lsp.py")
lsp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lsp)

# Deliberately small independent server implementation. Fragment every frame to
# exercise byte framing, not readline-based assumptions about JSON bodies.
SERVER = r'''
import json, os, pathlib, signal, subprocess, sys, time
ROOT = pathlib.Path(__file__).parent
RELOADED = False
(ROOT / "pid").write_text(str(os.getpid()))
def log(value):
    with (ROOT / "events.jsonl").open("a") as f:
        f.write(json.dumps(value) + "\n")
def send(value):
    body = json.dumps(dict(jsonrpc="2.0", **value), ensure_ascii=False).encode()
    data = b"Content-Type: application/vscode-jsonrpc; charset=utf-8\r\nContent-Length: " + str(len(body)).encode() + b"\r\n\r\n" + body
    for start in range(0, len(data), 7):
        sys.stdout.buffer.write(data[start:start + 7])
        sys.stdout.buffer.flush()
def recv():
    headers = {}
    while True:
        line = sys.stdin.buffer.readline()
        if not line: sys.exit(0)
        if line == b"\r\n": break
        key, value = line.split(b":", 1)
        headers[key.lower()] = value.strip()
    value = json.loads(sys.stdin.buffer.read(int(headers[b"content-length"])))
    log(value)
    return value

def response(request, result): send({"id": request["id"], "result": result})

request = recv()
assert request["method"] == "initialize"
log({"offline": os.environ.get("CARGO_NET_OFFLINE"), "auto_install": os.environ.get("RUSTUP_AUTO_INSTALL")})
if MODE == "malformed":
    sys.stdout.buffer.write(b"Content-Length: no\r\n\r\n")
    sys.stdout.buffer.flush()
    time.sleep(60)
if MODE == "initialize_error":
    send({"id": request["id"], "error": {"code": -32002, "message": "initialization failed"}})
else:
    response(request, {"capabilities": {"positionEncoding": "utf-8" if MODE == "encoding" else "utf-16",
                                       "definitionProvider": MODE != "unsupported",
                                       "referencesProvider": True, "implementationProvider": True}})
if MODE == "stall_write":
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    time.sleep(60)
while True:
    request = recv()
    method = request.get("method")
    if method == "exit": sys.exit(0)
    if method == "shutdown":
        if MODE in ("timeout", "stall_loading"): time.sleep(60)
        response(request, None)
    elif method == "textDocument/didOpen":
        uri = request["params"]["textDocument"]["uri"]
        if MODE == "exit_early": sys.exit(4)
        if MODE == "stall_loading": continue
        if MODE == "idle_only":
            send({"method": "experimental/serverStatus", "params": {"quiescent": True, "health": "ok"}})
            continue
        send({"method": "experimental/serverStatus", "params": {"quiescent": False, "health": "ok"}})
        requests = [
            ("workspace/configuration", {"items": [{"section": "rust-analyzer"}, {"section": "rust-analyzer.cargo.buildScripts.enable"}, {"section": "rust-analyzer.procMacro.enable"}, {"section": "unknown"}]}),
            ("client/registerCapability", {"registrations": []}),
            ("client/unregisterCapability", {"unregisterations": []}),
            ("window/workDoneProgress/create", {"token": "loading"}),
            ("workspace/workspaceFolders", {}),
            ("workspace/applyEdit", {"edit": {}}),
        ]
        for i, (method, params) in enumerate(requests):
            send({"id": "server-" + str(i), "method": method, "params": params})
            reply = recv()
            assert reply["id"] == "server-" + str(i), reply
        send({"method": "$/progress", "params": {"token": "loading", "value": {"kind": "end"}}})
        send({"method": "experimental/serverStatus", "params": {
            "quiescent": True, "health": "error" if MODE == "load_error" else "warning" if MODE == "warning" else "ok",
            "message": "offline dependency missing" if MODE in ("load_error", "warning") else None}})
    elif method in ("textDocument/definition", "textDocument/references", "textDocument/implementation"):
        if MODE == "reload" and not RELOADED:
            RELOADED = True
            send({"method": "experimental/serverStatus", "params": {"quiescent": False, "health": "ok"}})
            response(request, [])
            send({"method": "experimental/serverStatus", "params": {"quiescent": True, "health": "ok"}})
            continue
        if MODE == "timeout":
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
            time.sleep(60)
        elif MODE == "query_error":
            send({"id": request["id"], "error": {"code": -32801, "message": "content modified"}})
        elif MODE == "empty": response(request, [])
        elif MODE == "null": response(request, None)
        else:
            selection = {"start": {"line": 0, "character": 3}, "end": {"line": 0, "character": 7}}
            whole = {"start": {"line": 0, "character": 0}, "end": {"line": 0, "character": 7}}
            if MODE == "bad_position": selection["start"]["character"] = 1
            if MODE == "link":
                response(request, [{"targetUri": uri, "targetRange": whole, "targetSelectionRange": selection,
                                    "originSelectionRange": selection}])
            elif MODE == "external":
                response(request, [{"uri": "file:///outside/does%20not%20exist.rs", "range": selection}])
            elif MODE in ("ignored_result", "symlink_result"):
                target = (ROOT / ("ignored.rs" if MODE == "ignored_result" else "alias.rs")).as_uri()
                response(request, [{"uri": target, "range": selection}])
            else:
                value = {"uri": uri, "range": selection}
                response(request, value if MODE == "single" else [value])
'''


@unittest.skipUnless(os.name == "posix", "executable Python fixture uses POSIX shebangs")
class QueryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        subprocess.run(["git", "init", "-q", str(self.root)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.source = self.root / "source space.rs"
        self.source.write_text("😀 αname\nfn other() {}\n", encoding="utf-8")
        self.executable = self.root / "fake-lsp"

    def run_query(self, mode="location", **kwargs):
        self.executable.write_text("#!" + sys.executable + "\nMODE = " + repr(mode) + "\n" + SERVER)
        self.executable.chmod(0o755)
        arguments = dict(root=self.root, method="definition", path=self.source.name, line=1, column=3,
                         executable=str(self.executable), timeout=3)
        arguments.update(kwargs)
        result = lsp.query(**arguments)
        self.assert_dead()
        # The public result is always JSON serializable.
        json.dumps(result)
        return result

    def assert_dead(self):
        pid_file = self.root / "pid"
        if pid_file.exists():
            with self.assertRaises(ProcessLookupError):
                os.kill(int(pid_file.read_text()), 0)

    def events(self):
        return [json.loads(line) for line in (self.root / "events.jsonl").read_text().splitlines()]

    def test_location_and_unicode_codepoint_round_trip(self):
        result = self.run_query()
        self.assertEqual(result["status"], "ok", result)
        self.assertEqual(result["resolution_status"], "resolved")
        location = result["results"][0]
        self.assertEqual((location["line"], location["column"], location["end_column"]), (1, 3, 7))
        self.assertEqual(location["path"], self.source.name)
        self.assertEqual(location["uri"], self.source.as_uri())
        self.assertEqual(location["position_encoding"], "codepoint")
        query = next(e for e in self.events() if e.get("method") == "textDocument/definition")
        self.assertEqual(query["params"]["position"], {"line": 0, "character": 3})

    def test_location_link_preserves_all_ranges(self):
        location = self.run_query("link")["results"][0]
        self.assertEqual(location["range"]["start"], {"line": 1, "column": 1})
        self.assertEqual(location["target_selection_range"]["start"], {"line": 1, "column": 3})
        self.assertEqual(location["origin_selection_range"]["end"], {"line": 1, "column": 7})
        self.assertEqual(location["lsp_range"]["end"]["character"], 7)
        self.assertEqual(location["lsp_target_selection_range"]["start"]["character"], 3)

    def test_single_location(self):
        self.assertEqual(len(self.run_query("single")["results"]), 1)

    def test_references_and_implementation(self):
        for method in ("references", "implementation", "implementations"):
            with self.subTest(method=method):
                result = self.run_query(method=method)
                self.assertEqual(result["status"], "ok", result)
        request = next(e for e in self.events() if e.get("method") == "textDocument/references")
        self.assertEqual(request["params"]["context"], {"includeDeclaration": True})

    def test_startup_configuration_and_server_requests(self):
        self.assertEqual(self.run_query()["status"], "ok")
        events = self.events()
        initialize = events[0]["params"]
        config = initialize["initializationOptions"]
        self.assertFalse(config["check"]["enable"])
        self.assertFalse(config["checkOnSave"])
        self.assertFalse(config["cargo"]["buildScripts"]["enable"])
        self.assertFalse(config["procMacro"]["enable"])
        self.assertTrue(config["cargo"]["offline"])
        self.assertEqual(config["cargo"]["extraEnv"]["CARGO_NET_OFFLINE"], "true")
        self.assertIsNone(config["cargo"]["sysroot"])
        self.assertTrue(initialize["capabilities"]["experimental"]["serverStatusNotification"])
        self.assertIn({"offline": "true", "auto_install": "0"}, events)
        replies = {e.get("id"): e for e in events if "method" not in e}
        self.assertEqual(replies["server-0"]["result"], [config, False, False, None])
        self.assertIsNone(replies["server-1"]["result"])
        self.assertIsNone(replies["server-2"]["result"])
        self.assertIsNone(replies["server-3"]["result"])
        self.assertEqual(replies["server-4"]["result"][0]["uri"], self.root.as_uri())
        self.assertEqual(replies["server-5"]["error"]["code"], -32601)
        methods = [e.get("method") for e in events]
        self.assertLess(methods.index("initialized"), methods.index("textDocument/didOpen"))
        self.assertLess(methods.index("textDocument/didOpen"), methods.index("textDocument/definition"))
        self.assertIn("shutdown", methods)
        self.assertIn("exit", methods)

    def test_empty_only_after_server_ready(self):
        for mode in ("empty", "null"):
            result = self.run_query(mode)
            self.assertEqual(result["status"], "ok", result)
            self.assertEqual(result["resolution_status"], "no_results")
            self.assertEqual(result["results"], [])

    def test_reload_racing_empty_response_is_retried(self):
        result = self.run_query("reload")
        self.assertEqual(result["resolution_status"], "resolved", result)
        queries = [e for e in self.events() if e.get("method") == "textDocument/definition"]
        self.assertEqual(len(queries), 2)

    def test_idle_without_observed_loading_is_not_ready(self):
        result = self.run_query("idle_only", timeout=0.4)
        self.assertEqual(result["error"]["code"], "timeout", result)
        self.assertEqual(result["error"]["phase"], "project_loading")
        self.assertFalse(any(e.get("method") == "textDocument/definition" for e in self.events()))

    def test_external_uri_retained_without_fabricated_codepoints(self):
        location = self.run_query("external")["results"][0]
        self.assertIsNone(location["path"])
        self.assertEqual(location["uri"], "file:///outside/does%20not%20exist.rs")
        self.assertEqual(location["position_encoding"], "utf-16")
        self.assertEqual(location["column"], 4)

    def test_ignored_and_symlink_results_not_read(self):
        (self.root / ".gitignore").write_text("ignored.rs\n")
        (self.root / "ignored.rs").write_text("not long enough")
        (self.root / "alias.rs").symlink_to(self.source)
        for mode in ("ignored_result", "symlink_result"):
            location = self.run_query(mode)["results"][0]
            self.assertEqual(location["position_encoding"], "utf-16")
            self.assertEqual(location["column"], 4)

    def test_explicit_errors(self):
        for mode, code in [("initialize_error", "server_error"), ("query_error", "server_error"),
                           ("load_error", "project_load_error"), ("encoding", "unsupported_encoding"),
                           ("unsupported", "unsupported_method"), ("malformed", "protocol_error"),
                           ("exit_early", "server_exited"), ("bad_position", "protocol_error")]:
            with self.subTest(mode=mode):
                result = self.run_query(mode)
                self.assertEqual(result["status"], "error", result)
                self.assertEqual(result["error"]["code"], code, result)
                self.assertEqual(result["results"], [])

    def test_warning_is_preserved_with_limitations(self):
        result = self.run_query("warning")
        self.assertEqual(result["status"], "ok", result)
        self.assertEqual(result["server_health"], "warning")
        self.assertEqual(result["server_message"], "offline dependency missing")
        self.assertTrue(result["limitations"])

    def test_timeouts_are_bounded_and_process_is_reaped(self):
        for mode, phase in (("stall_loading", "project_loading"), ("timeout", "query")):
            with self.subTest(mode=mode):
                start = time.monotonic()
                result = self.run_query(mode, timeout=0.6)
                self.assertLess(time.monotonic() - start, 2)
                self.assertEqual(result["error"]["code"], "timeout", result)
                self.assertEqual(result["error"]["phase"], phase)

    def test_large_write_timeout_is_bounded(self):
        self.source.write_text("a" * 2000000)
        start = time.monotonic()
        result = self.run_query("stall_write", timeout=0.3)
        self.assertLess(time.monotonic() - start, 2)
        self.assertEqual(result["error"]["code"], "timeout", result)

    def test_missing_executable_is_not_lexical_fallback(self):
        result = lsp.query(self.root, "definition", self.source.name, 1, 1,
                           executable=str(self.root / "missing"), timeout=1)
        self.assertEqual(result["error"]["code"], "unavailable")
        self.assertEqual(result["mode"], "semantic")
        self.assertEqual(result["results"], [])

    def test_invalid_inputs_never_start_server(self):
        (self.root / ".gitignore").write_text("ignored.rs\n")
        (self.root / "ignored.rs").write_text("fn ignored() {}")
        (self.root / "alias.rs").symlink_to(self.source)
        (self.root / "source.txt").write_text("fn invalid() {}")
        (self.root / "linked").symlink_to(self.root, target_is_directory=True)
        invalid = [{"path": "../outside.rs"}, {"path": str(self.source)}, {"path": "ignored.rs"},
                   {"path": "alias.rs"}, {"path": "linked/" + self.source.name}, {"path": "source.txt"},
                   {"line": 0}, {"column": 20}, {"column": True}, {"method": "hover"},
                   {"timeout": float("inf")}, {"timeout": 0}]
        for arguments in invalid:
            with self.subTest(arguments=arguments):
                result = self.run_query(**arguments)
                self.assertEqual(result["status"], "error", result)
                self.assertFalse((self.root / "pid").exists())

    def test_non_git_root_is_explicit_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "lib.rs").write_text("fn main() {}")
            result = lsp.query(root, "definition", "lib.rs", 1, 1)
            self.assertEqual(result["error"]["code"], "invalid_root")

    def test_line_breaks_and_surrogate_conversion(self):
        self.assertEqual(lsp._lines("a\u2028b\r\nc\n"), ["a\u2028b", "c", ""])
        self.assertEqual(lsp._point({"line": 0, "character": 2}, ["😀x"]), {"line": 1, "column": 2})
        with self.assertRaises(lsp._Failure):
            lsp._point({"line": 0, "character": 1}, ["😀x"])


if __name__ == "__main__":
    unittest.main()
