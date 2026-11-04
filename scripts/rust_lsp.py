#!/usr/bin/env python3
"""Optional, bounded rust-analyzer LSP queries; no lexical fallback or installation.

Positions supplied to query() are 1-based Unicode codepoints. Results preserve raw
LSP ranges as well as converted 1-based ranges. External/unreadable source is not
opened: those locations are explicitly labelled UTF-16 instead of guessing.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import queue
import signal
import subprocess
import threading
import time
from urllib.parse import unquote, urlsplit

_METHODS = {"definition": "textDocument/definition", "references": "textDocument/references",
            "implementation": "textDocument/implementation", "implementations": "textDocument/implementation"}
_EXCLUDED = {"target", "build", "dist", "node_modules", ".git", ".venv"}
_CONFIG = {
    "checkOnSave": False,
    "check": {"enable": False},
    "cargo": {"offline": True, "buildScripts": {"enable": False},
              "extraEnv": {"CARGO_NET_OFFLINE": "true"}, "sysroot": None},
    "procMacro": {"enable": False},
}
_LIMITATIONS = [
    "No lexical fallback; semantic results depend on rust-analyzer and available offline dependencies.",
    "Checks, build scripts and procedural macros are disabled; generated/cfg-dependent code may be incomplete.",
    "Sysroot discovery is disabled; standard-library resolution may be unavailable.",
    "External, ignored and symlinked sources are not read; their positions remain explicitly UTF-16.",
    "A ready server and an empty result do not prove that no definition/reference/implementation exists.",
]


class _Failure(Exception):
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


def _remaining(deadline):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise _Failure("timeout", "rust-analyzer deadline exceeded")
    return remaining


class _Client:
    """One reader and writer avoid blocking indefinitely on either subprocess pipe."""
    def __init__(self, executable, root, deadline):
        self.deadline = deadline
        self.incoming = queue.Queue(maxsize=1024)
        self.outgoing = queue.Queue()
        self.counter = 0
        self.ready = False
        self.saw_loading = False
        self.health = None
        self.status_message = None
        self.stderr = bytearray()
        env = dict(os.environ, CARGO_NET_OFFLINE="true", RUSTUP_AUTO_INSTALL="0",
                   RUSTUP_NO_UPDATE_CHECK="1")
        self.process = subprocess.Popen(
            [os.fspath(executable)], cwd=root, env=env, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=(os.name == "posix"),
        )
        self.threads = [threading.Thread(target=target, daemon=True)
                        for target in (self._read, self._write, self._stderr)]
        for thread in self.threads:
            thread.start()

    def _read(self):
        try:
            stream = self.process.stdout
            while True:
                headers, size = {}, 0
                while True:
                    line = stream.readline(8193)
                    size += len(line)
                    if not line:
                        raise _Failure("server_exited", "rust-analyzer closed stdout")
                    if size > 32768 or len(line) > 8192:
                        raise _Failure("protocol_error", "oversized LSP header")
                    if line in (b"\r\n", b"\n"):
                        break
                    key, sep, value = line.partition(b":")
                    if not sep or key.lower() in headers:
                        raise _Failure("protocol_error", "malformed or duplicate LSP header")
                    headers[key.lower()] = value.strip()
                try:
                    length = int(headers[b"content-length"])
                except (KeyError, ValueError):
                    raise _Failure("protocol_error", "invalid Content-Length") from None
                if not 0 < length <= 16 * 1024 * 1024:
                    raise _Failure("protocol_error", "LSP message size out of bounds")
                chunks, remaining = [], length
                while remaining:
                    chunk = stream.read(remaining)
                    if not chunk:
                        raise _Failure("protocol_error", "truncated LSP body")
                    chunks.append(chunk)
                    remaining -= len(chunk)
                message = json.loads(b"".join(chunks))
                if not isinstance(message, dict) or message.get("jsonrpc") != "2.0":
                    raise _Failure("protocol_error", "invalid JSON-RPC envelope")
                self.incoming.put(message)
        except Exception as exc:
            self.incoming.put(exc if isinstance(exc, _Failure) else _Failure("protocol_error", str(exc)))

    def _write(self):
        while True:
            item = self.outgoing.get()
            if item is None:
                return
            message, done, errors = item
            try:
                body = json.dumps(message, ensure_ascii=False).encode("utf-8")
                self.process.stdin.write(f"Content-Length: {len(body)}\r\n\r\n".encode("ascii") + body)
                self.process.stdin.flush()
            except Exception as exc:
                errors.append(exc)
            finally:
                done.set()

    def _stderr(self):
        while True:
            chunk = self.process.stderr.read(4096)
            if not chunk:
                return
            self.stderr.extend(chunk)
            del self.stderr[:-8192]

    def send(self, message):
        done, errors = threading.Event(), []
        self.outgoing.put((dict(jsonrpc="2.0", **message), done, errors))
        if not done.wait(_remaining(self.deadline)):
            raise _Failure("timeout", "timeout writing to rust-analyzer")
        if errors:
            raise _Failure("server_exited", f"cannot write to rust-analyzer: {errors[0]}")

    def notify(self, method, params):
        self.send({"method": method, "params": params})

    def _configuration(self, item):
        section = item.get("section", "")
        if section in ("", "rust-analyzer"):
            return _CONFIG
        if section.startswith("rust-analyzer."):
            section = section[len("rust-analyzer."):]
        value = _CONFIG
        for part in section.split("."):
            if not isinstance(value, dict):
                return None
            value = value.get(part)
        return value

    def next(self):
        try:
            message = self.incoming.get(timeout=_remaining(self.deadline))
        except queue.Empty:
            raise _Failure("timeout", "timeout waiting for rust-analyzer") from None
        if isinstance(message, Exception):
            raise message
        method = message.get("method")
        if method and "id" in message:
            if method == "workspace/configuration":
                result = [self._configuration(item) for item in message.get("params", {}).get("items", [])]
            elif method in {"client/registerCapability", "client/unregisterCapability", "window/workDoneProgress/create"}:
                result = None
            elif method == "workspace/workspaceFolders":
                result = self.folders
            else:
                self.send({"id": message["id"], "error": {"code": -32601, "message": "Unsupported client request"}})
                return {}
            self.send({"id": message["id"], "result": result})
            return {}
        if method == "experimental/serverStatus":
            params = message.get("params", {})
            if params.get("quiescent") is False:
                self.saw_loading = True
            self.ready = self.saw_loading and params.get("quiescent") is True
            self.health = params.get("health")
            self.status_message = params.get("message")
        return message

    def request(self, method, params):
        self.counter += 1
        request_id = self.counter
        self.send({"id": request_id, "method": method, "params": params})
        while True:
            message = self.next()
            if message.get("id") == request_id and "method" not in message:
                if "error" in message:
                    raise _Failure("server_error", json.dumps(message["error"], ensure_ascii=False))
                if "result" not in message:
                    raise _Failure("protocol_error", "response has neither result nor error")
                return message["result"]

    def close(self):
        # Cleanup has a separate small budget even if the query exhausted its deadline.
        self.deadline = time.monotonic() + 0.25
        if self.process.poll() is None:
            try:
                self.request("shutdown", None)
                self.notify("exit", None)
                self.process.wait(timeout=0.15)
            except Exception:
                pass
        # Kill the dedicated process group too: a server can leave pipe-holding children.
        if os.name == "posix":
            try:
                os.killpg(self.process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        elif self.process.poll() is None:
            self.process.terminate()
        try:
            self.process.wait(timeout=0.2)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait(timeout=0.2)
        if os.name == "posix":
            try:
                os.killpg(self.process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        self.outgoing.put(None)
        for thread in self.threads:
            thread.join(timeout=0.1)
        for stream in (self.process.stdin, self.process.stdout, self.process.stderr):
            try:
                stream.close()
            except OSError:
                # Buffered stdin can flush into a server killed after a write timeout.
                pass


def _visible_source(root, relative, deadline):
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts or candidate.suffix != ".rs":
        raise _Failure("invalid_path", "path must be a repository-relative .rs source without traversal")
    if _EXCLUDED.intersection(candidate.parts):
        raise _Failure("invalid_path", "excluded source path")
    current = root
    for part in candidate.parts:
        current = current / part
        if current.is_symlink():
            raise _Failure("invalid_path", "symlinked source paths are not allowed")
    if not current.is_file() or not current.resolve().is_relative_to(root):
        raise _Failure("invalid_path", "source is missing or outside the repository")
    try:
        result = subprocess.run(
            ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", ":(literal)" + candidate.as_posix()],
            cwd=root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=_remaining(deadline), check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise _Failure("invalid_root", "source validation requires a Git worktree") from exc
    if os.fsencode(candidate.as_posix()) not in result.stdout.split(b"\0"):
        raise _Failure("invalid_path", "ignored/non-source path is not Git-visible")
    return current


def _lines(text):
    # splitlines() incorrectly treats Unicode separators as LSP line breaks.
    return text.replace("\r\n", "\n").split("\n")


def _point(point, lines):
    line, units = point["line"], point["character"]
    if not isinstance(line, int) or not isinstance(units, int) or line < 0 or units < 0:
        raise _Failure("protocol_error", "invalid LSP position")
    if lines is None:
        return {"line": line + 1, "column": units + 1}
    if line >= len(lines):
        raise _Failure("protocol_error", "LSP position exceeds source lines")
    count = 0
    for index, char in enumerate(lines[line]):
        if count == units:
            return {"line": line + 1, "column": index + 1}
        count += 2 if ord(char) > 0xFFFF else 1
        if count > units:
            raise _Failure("protocol_error", "LSP position splits a UTF-16 surrogate pair")
    if count != units:
        raise _Failure("protocol_error", "LSP position exceeds source line")
    return {"line": line + 1, "column": len(lines[line]) + 1}


def _range(value, lines):
    return {"start": _point(value["start"], lines), "end": _point(value["end"], lines)}


def _location(value, root, source_path, source_lines, deadline):
    _remaining(deadline)
    link = "targetUri" in value
    uri = value["targetUri"] if link else value["uri"]
    raw_range = value["targetRange"] if link else value["range"]
    selection = value.get("targetSelectionRange", raw_range) if link else raw_range
    parsed = urlsplit(uri)
    relative, lines = None, None
    if parsed.scheme == "file" and parsed.netloc in ("", "localhost"):
        local = Path(unquote(parsed.path))
        if local.is_absolute() and ".." not in local.parts:
            try:
                relative = local.relative_to(root).as_posix()
            except ValueError:
                pass
            if relative is not None:
                try:
                    safe = _visible_source(root, relative, deadline)
                    lines = source_lines if safe == source_path else _lines(safe.read_text(encoding="utf-8"))
                except _Failure as exc:
                    if exc.code == "timeout":
                        raise
                except (OSError, UnicodeError):
                    pass
    converted = _range(selection, lines)
    result = {"uri": uri, "path": relative, "position_encoding": "codepoint" if lines is not None else "utf-16",
              "range": _range(raw_range, lines), "lsp_range": raw_range,
              **converted["start"], "end_line": converted["end"]["line"], "end_column": converted["end"]["column"]}
    if link:
        result.update(target_selection_range=converted, lsp_target_selection_range=selection)
        if "originSelectionRange" in value:
            result["origin_selection_range"] = _range(value["originSelectionRange"], source_lines)
            result["lsp_origin_selection_range"] = value["originSelectionRange"]
    return result


def query(root: Path, method: str, path: str, line: int, column: int,
          executable="rust-analyzer", timeout=30) -> dict:
    """Run one optional semantic query. All failures are explicit JSON-compatible data.

    Requires an observed experimental/serverStatus loading-to-quiescent cycle,
    not a sleep or an empty probe. Servers lacking it fail explicitly on timeout.
    The timeout covers validation, startup/loading and request, plus bounded cleanup.
    Only a preinstalled executable is run; no Cargo build/check/install is invoked.
    """
    result = {"schema_version": 1, "mode": "semantic", "semantics": "rust-analyzer",
              "status": "error", "resolution_status": "error", "results": [],
              "limitations": list(_LIMITATIONS), "analyzer": os.fspath(executable), "method": method}
    client = None
    phase = "validation"
    try:
        if not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or timeout <= 0:
            raise _Failure("invalid_argument", "timeout must be finite and positive")
        deadline = time.monotonic() + timeout
        root = Path(root).resolve(strict=True)
        if method not in _METHODS:
            raise _Failure("invalid_argument", "method must be definition, references or implementation")
        source = _visible_source(root, path, deadline)
        text = source.read_text(encoding="utf-8")
        lines = _lines(text)
        if (type(line) is not int or type(column) is not int or line < 1 or line > len(lines)
                or column < 1 or column > len(lines[line - 1]) + 1):
            raise _Failure("invalid_position", "line/column must be 1-based codepoints within source")
        position = {"line": line - 1, "character": len(lines[line - 1][:column - 1].encode("utf-16-le")) // 2}
        phase = "initialize"
        client = _Client(executable, root, deadline)
        client.folders = [{"uri": root.as_uri(), "name": root.name}]
        initialized = client.request("initialize", {
            "processId": os.getpid(), "rootUri": root.as_uri(), "workspaceFolders": client.folders,
            "capabilities": {"general": {"positionEncodings": ["utf-16"]},
                             "workspace": {"configuration": True, "workspaceFolders": True},
                             "window": {"workDoneProgress": True},
                             "textDocument": {"definition": {"linkSupport": True}, "implementation": {"linkSupport": True}},
                             "experimental": {"serverStatusNotification": True}},
            "initializationOptions": _CONFIG,
        })
        capabilities = initialized.get("capabilities", {})
        if capabilities.get("positionEncoding", "utf-16") != "utf-16":
            raise _Failure("unsupported_encoding", "server did not select UTF-16")
        provider = {"definition": "definitionProvider", "references": "referencesProvider",
                    "implementation": "implementationProvider", "implementations": "implementationProvider"}[method]
        if capabilities.get(provider) is False:
            raise _Failure("unsupported_method", f"server does not support {method}")
        client.notify("initialized", {})
        client.notify("workspace/didChangeConfiguration", {"settings": {"rust-analyzer": _CONFIG}})
        client.notify("textDocument/didOpen", {"textDocument": {"uri": source.as_uri(), "languageId": "rust", "version": 1, "text": text}})
        phase = "project_loading"
        while not client.ready:
            client.next()
        if client.health == "error":
            raise _Failure("project_load_error", client.status_message or "rust-analyzer project loading failed")
        result["server_health"] = client.health
        result["server_message"] = client.status_message
        phase = "query"
        params = {"textDocument": {"uri": source.as_uri()}, "position": position}
        if method == "references":
            params["context"] = {"includeDeclaration": True}
        while True:
            values = client.request(_METHODS[method], params)
            if client.ready:
                break
            # A request can race a reload. Never accept its startup-empty reply.
            phase = "project_loading"
            while not client.ready:
                client.next()
            if client.health == "error":
                raise _Failure("project_load_error", client.status_message or "project reload failed")
            phase = "query"
        result["server_health"] = client.health
        result["server_message"] = client.status_message
        values = [] if values is None else values if isinstance(values, list) else [values]
        result["results"] = [_location(value, root, source, lines, deadline) for value in values]
        result.update(status="ok", resolution_status="resolved" if values else "no_results")
    except Exception as exc:
        code = exc.code if isinstance(exc, _Failure) else "timeout" if isinstance(exc, subprocess.TimeoutExpired) else "adapter_error"
        if isinstance(exc, FileNotFoundError):
            code = "unavailable" if phase == "initialize" else "invalid_path"
        result.update(status="error", resolution_status="error", results=[],
                      error={"code": code, "phase": phase, "message": str(exc)})
    finally:
        if client is not None:
            client.close()
    return result
