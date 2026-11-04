#!/usr/bin/env python3
"""Reproducible repair tasks in disposable source snapshots, never the live checkout."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fnmatch
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
SUITE = ROOT / ".agents/evals/tasks.json"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def isolated_environment():
    # Repository-selection/config overrides must never redirect snapshot Git
    # commands (including commands run by acceptance tests) to the live checkout.
    return {key: value for key, value in os.environ.items() if not key.startswith("GIT_") and key not in {"PYTHONPATH", "PYTHONHOME", "CARGO_TARGET_DIR"}}


def git(root, *args):
    return subprocess.check_output(["git", "-c", "core.hooksPath=/dev/null", "-c", "commit.gpgsign=false", *args], cwd=root, env=isolated_environment())


def tasks():
    suite = json.loads(SUITE.read_text())
    if suite["schema_version"] != 1:
        raise ValueError("Unsupported evaluation suite schema")
    result = {}
    for task in suite["tasks"]:
        if task["id"] in result or not task["acceptance"] or not task["expected_failure_contains"]:
            raise ValueError("Duplicate or incomplete evaluation task")
        if not 1 <= task["timeout_seconds"] <= 600:
            raise ValueError("Acceptance timeout must be 1..600 seconds")
        result[task["id"]] = task
    return result


def run_directory(run_id):
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", run_id):
        raise ValueError("Run id must contain 1..64 lowercase letters, digits or hyphens")
    expected = ROOT / "outputs/agent-evals" / run_id
    path = expected.resolve()
    if path != expected.absolute():
        raise ValueError("Evaluation directories must not contain symlink components")
    return path


def local_path(workspace, name):
    path = (workspace / name).resolve()
    if not path.is_relative_to(workspace.resolve()) or Path(name).is_absolute():
        raise ValueError(f"Path outside evaluation workspace: {name}")
    return path


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def load_run(run_id):
    directory = run_directory(run_id)
    if (directory / "workspace").is_symlink():
        raise ValueError("Evaluation workspace must not be a symlink")
    state = json.loads((directory / "run.json").read_text())
    if state["run_id"] != run_id or state["task_sha256"] != digest(state["task"]):
        raise ValueError("Run metadata does not match its task definition")
    return directory, state


def save_run(directory, state):
    write_json(directory / "run.json", state)


def changed_paths(workspace, baseline="HEAD"):
    result = set()
    for arguments in (("diff", "--no-renames", "--name-only", "-z", baseline, "--"), ("ls-files", "--others", "--exclude-standard", "-z")):
        result.update(os.fsdecode(path) for path in git(workspace, *arguments).split(b"\0") if path)
    return sorted(result)


def prepare(task_id, run_id, variant, ref="HEAD", agent_metadata=None):
    task = tasks()[task_id]
    directory = run_directory(run_id)
    if directory.exists():
        raise ValueError("Run directory already exists; choose a fresh run id")
    commit = git(ROOT, "rev-parse", "--verify", "--end-of-options", ref + "^{commit}").decode().strip()
    archive = git(ROOT, "archive", commit)
    # Validate every member before writing; source snapshots must not contain links.
    with tarfile.open(fileobj=io.BytesIO(archive)) as source:
        for member in source.getmembers():
            name = Path(member.name)
            if name.is_absolute() or ".." in name.parts or not (member.isfile() or member.isdir()):
                raise ValueError("Unsafe or unsupported Git archive member")
        directory.mkdir(parents=True)
        workspace = directory / "workspace"
        workspace.mkdir()
        source.extractall(workspace, filter="data")
    for mutation in task["mutations"]:
        path = local_path(workspace, mutation["path"])
        content = path.read_text()
        if content.count(mutation["old"]) != 1:
            raise ValueError(f"Task seed does not match {mutation['path']} at {commit}; no live source changed")
        path.write_text(content.replace(mutation["old"], mutation["new"], 1))
    git(workspace, "init", "--quiet")
    git(workspace, "add", ".")
    git(workspace, "-c", "user.name=Agent evaluation", "-c", "user.email=eval@example.invalid", "commit", "--quiet", "-m", "Seeded repair task")
    state = {"schema_version": 1, "run_id": run_id, "task": task, "task_sha256": digest(task), "baseline_commit": commit, "seed_commit": git(workspace, "rev-parse", "HEAD").decode().strip(),
             "variant": variant, "agent_metadata": agent_metadata or {"agent": "unspecified", "model": "unspecified", "settings": {}}, "created_at": datetime.now(timezone.utc).isoformat(), "seed_verified": False,
             "started_unix": None, "instrumented_context_bytes": 0, "retrieval_events": []}
    save_run(directory, state)
    (directory / "TASK.md").write_text(f"# {task['id']}\n\n{task['brief']}\n\nAllowed source edits: {', '.join(task['allowed_paths'])}.\nDo not change acceptance tests or the evaluator.\n\nThis is a controlled repair task, not a security sandbox or a blind benchmark.\n")
    return {"run_id": run_id, "workspace": str(workspace), "brief": task["brief"], "tier": task["tier"],
            "next": f"Run verify {run_id} before handing the workspace to an agent; timing starts only after the expected seeded failure is observed."}


def execute_checks(directory, state, prefix):
    workspace = directory / "workspace"
    results, outputs = [], []
    for index, command in enumerate(state["task"]["acceptance"]):
        argv = [sys.executable if arg == "{python}" else arg for arg in command]
        log = directory / f"{prefix}-{index}.log"
        started = time.perf_counter()
        timed_out = False
        with log.open("x") as stream:
            try:
                process = subprocess.Popen(argv, cwd=workspace, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True, env=dict(isolated_environment(), CARGO_NET_OFFLINE="true", PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(directory / f"{prefix}-{index}-pycache")))
                try:
                    returncode = process.wait(timeout=state["task"]["timeout_seconds"])
                except subprocess.TimeoutExpired:
                    timed_out = True
                    os.killpg(process.pid, signal.SIGKILL)
                    returncode = process.wait()
            except OSError as exc:
                stream.write(str(exc))
                returncode = 127
        outputs.append(log.read_text(errors="replace"))
        results.append({"argv": argv, "returncode": returncode, "timed_out": timed_out, "seconds": time.perf_counter() - started, "log": log.name})
    return results, "\n".join(outputs)


def verify(run_id):
    directory, state = load_run(run_id)
    if state["seed_verified"]:
        raise ValueError("Seed already verified; use a new run for a fresh timed attempt")
    if changed_paths(directory / "workspace", state["seed_commit"]):
        raise ValueError("Verify the seed before editing the evaluation workspace")
    results, output = execute_checks(directory, state, "seed")
    markers = state["task"]["expected_failure_contains"]
    valid = any(item["returncode"] != 0 and not item["timed_out"] for item in results) and all(marker in output for marker in markers)
    state["seed_checks"] = results
    state["seed_verified"] = valid
    state["started_unix"] = time.time() if valid else None
    save_run(directory, state)
    return {"run_id": run_id, "seed_verified": valid, "checks": results,
            "note": "A missing compiler/dependency or unrelated failure does not establish a valid repair task."}


def record_context(run_id, text, source):
    directory, state = load_run(run_id)
    if not state["seed_verified"]:
        raise ValueError("Verify the seed before measuring agent retrieval")
    size = len(text.encode())
    state["instrumented_context_bytes"] += size
    state["retrieval_events"].append({"source": source, "bytes": size})
    save_run(directory, state)
    return text


def read_context(run_id, path, offset=1, limit=100):
    if offset < 1 or not 1 <= limit <= 500:
        raise ValueError("Use a positive offset and a limit of 1..500 lines")
    directory, _ = load_run(run_id)
    workspace = directory / "workspace"
    file = local_path(workspace, path)
    text = "\n".join(f"{index + offset}: {line}" for index, line in enumerate(file.read_text().splitlines()[offset - 1:offset - 1 + limit])) + "\n"
    return record_context(run_id, text, f"{path}:{offset}-{offset + limit - 1}")


def navigation_context(run_id, component, limit=25):
    if not 1 <= limit <= 100:
        raise ValueError("Context limit must be 1..100 files")
    directory, _ = load_run(run_id)
    workspace = directory / "workspace"
    result = subprocess.run([sys.executable, str(workspace / "scripts/agent_tools.py"), "context", component, "--limit", str(limit)], cwd=workspace, text=True, capture_output=True, timeout=30, env=isolated_environment())
    if result.returncode:
        raise ValueError(result.stderr.strip() or "Context retrieval failed")
    return record_context(run_id, result.stdout, "context:" + component)


def usage_record(path):
    if path is None:
        return {"input_tokens": None, "output_tokens": None, "source": "not provided"}
    usage = json.loads(Path(path).read_text())
    for key in ("input_tokens", "output_tokens"):
        if type(usage.get(key)) is not int or usage[key] < 0:
            raise ValueError("Token usage requires nonnegative integer input_tokens and output_tokens")
    if not isinstance(usage.get("source"), str) or not usage["source"].strip():
        raise ValueError("Provide the provenance of exported token usage in source")
    return {key: usage[key] for key in ("input_tokens", "output_tokens", "source")}


def grade(run_id, usage_path=None):
    directory, state = load_run(run_id)
    if not state["seed_verified"]:
        raise ValueError("A valid seeded failure must be verified before grading")
    usage = usage_record(usage_path)
    elapsed = time.time() - state["started_unix"]
    if not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError("Clock changed; repeat with a fresh run")
    changed = changed_paths(directory / "workspace", state["seed_commit"])
    violations = [path for path in changed if not any(fnmatch.fnmatchcase(path, pattern) for pattern in state["task"]["allowed_paths"])]
    attempt = len(list(directory.glob("report-*.json"))) + 1
    results, _ = execute_checks(directory, state, f"attempt-{attempt:03}")
    report = {"schema_version": 1, "run_id": run_id, "task_id": state["task"]["id"], "task_sha256": state["task_sha256"],
              "baseline_commit": state["baseline_commit"], "variant": state["variant"], "agent_metadata": state["agent_metadata"], "attempt": attempt,
              "passed": not violations and bool(changed) and all(check["returncode"] == 0 for check in results),
              "elapsed_seconds": elapsed, "instrumented_context_bytes": state["instrumented_context_bytes"], "usage": usage,
              "changed_paths": changed, "scope_violations": violations, "checks": results,
              "limits": "Wall time includes idle time; bytes cover only evaluator-mediated retrieval. Exported token counts are caller-supplied, never inferred. Scope/tests are a quality check, not a security sandbox."}
    write_json(directory / f"report-{attempt:03}.json", report)
    return report


def compare(paths):
    reports = [json.loads(Path(path).read_text()) for path in paths]
    pairs = {(r["baseline_commit"], r["task_sha256"], digest(r["agent_metadata"])) for r in reports}
    if len(pairs) != 1:
        raise ValueError("Compare the same baseline commit, task definition and agent/model/settings; do not pool unlike trials")
    if len({r["run_id"] for r in reports}) != len(reports):
        raise ValueError("Use one final attempt per independent run")
    groups = {}
    for report in reports:
        groups.setdefault(report["variant"], []).append(report)
    result = {}
    for variant, items in sorted(groups.items()):
        measured = [item for item in items if item["usage"]["input_tokens"] is not None]
        result[variant] = {"runs": len(items), "passed": sum(item["passed"] for item in items),
                           "median_elapsed_seconds": statistics.median(item["elapsed_seconds"] for item in items),
                           "median_instrumented_context_bytes": statistics.median(item["instrumented_context_bytes"] for item in items),
                           "runs_with_token_usage": len(measured),
                           "median_input_tokens": statistics.median(item["usage"]["input_tokens"] for item in measured) if measured else None,
                           "median_output_tokens": statistics.median(item["usage"]["output_tokens"] for item in measured) if measured else None}
    return {"agent_metadata": reports[0]["agent_metadata"], "variants": result, "note": "Descriptive observations only. Match agent/model/settings/budget and repeat paired trials before making efficiency claims."}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    new = sub.add_parser("prepare")
    new.add_argument("task")
    new.add_argument("run_id")
    new.add_argument("--variant", required=True)
    new.add_argument("--ref", default="HEAD")
    new.add_argument("--agent", default="unspecified", help="Agent client/version label, not an agent launcher")
    new.add_argument("--model", default="unspecified")
    new.add_argument("--settings", type=Path, help="JSON object with seed, reasoning effort, budget and other controlled settings")
    seed = sub.add_parser("verify")
    seed.add_argument("run_id")
    read = sub.add_parser("read")
    read.add_argument("run_id")
    read.add_argument("path")
    read.add_argument("--offset", type=int, default=1)
    read.add_argument("--limit", type=int, default=100)
    navigation = sub.add_parser("context")
    navigation.add_argument("run_id")
    navigation.add_argument("component")
    navigation.add_argument("--limit", type=int, default=25)
    check = sub.add_parser("grade")
    check.add_argument("run_id")
    check.add_argument("--usage", type=Path)
    summary = sub.add_parser("compare")
    summary.add_argument("reports", nargs="+", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            result = [{key: task[key] for key in ("id", "tier", "brief", "allowed_paths")} for task in tasks().values()]
        elif args.command == "prepare":
            settings = {} if args.settings is None else json.loads(args.settings.read_text())
            if not isinstance(settings, dict):
                raise ValueError("Agent settings must be a JSON object")
            result = prepare(args.task, args.run_id, args.variant, args.ref, {"agent": args.agent, "model": args.model, "settings": settings})
        elif args.command == "verify":
            result = verify(args.run_id)
        elif args.command == "read":
            print(read_context(args.run_id, args.path, args.offset, args.limit), end="")
            return 0
        elif args.command == "context":
            print(navigation_context(args.run_id, args.component, args.limit), end="")
            return 0
        elif args.command == "grade":
            result = grade(args.run_id, args.usage)
        else:
            result = compare(args.reports)
        print(json.dumps(result, indent=2, sort_keys=True))
        return 1 if (args.command == "grade" and not result["passed"]) or (args.command == "verify" and not result["seed_verified"]) else 0
    except (OSError, ValueError, KeyError, subprocess.SubprocessError, tarfile.TarError) as exc:
        print(f"agent_eval: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
