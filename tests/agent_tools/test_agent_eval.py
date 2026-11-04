"""Exercise evaluation lifecycle without launching an agent, network or native build."""
import contextlib
import copy
import importlib.util
import io
import json
import os
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("agent_eval", ROOT / "scripts/agent_eval.py")
tools = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tools)


class EvaluationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.suite = self.root / ".agents/evals/tasks.json"
        self.suite.parent.mkdir(parents=True)
        self.task = {
            "id": "tiny-repair", "tier": "offline", "brief": "Restore the documented value.",
            "allowed_paths": ["app.py"], "mutations": [{"path": "app.py", "old": "VALUE = 1", "new": "VALUE = 2"}],
            "acceptance": [["{python}", "-c", "from app import VALUE; assert VALUE == 1, 'EVAL_EXPECTED_FAILURE'"]],
            "expected_failure_contains": ["EVAL_EXPECTED_FAILURE"], "timeout_seconds": 5,
        }
        self.suite.write_text(json.dumps({"schema_version": 1, "tasks": [self.task]}))
        (self.root / "app.py").write_text("VALUE = 1\n")
        (self.root / ".gitignore").write_text("/outputs/\n__pycache__/\n")
        tools.git(self.root, "init", "--quiet")
        tools.git(self.root, "add", ".")
        tools.git(self.root, "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "--quiet", "-m", "fixture")
        self.addCleanup(patch.stopall)
        patch.object(tools, "ROOT", self.root).start()
        patch.object(tools, "SUITE", self.suite).start()

    def prepare(self, run="trial", variant="baseline"):
        return tools.prepare("tiny-repair", run, variant)

    def workspace(self, run="trial"):
        return tools.run_directory(run) / "workspace"

    def test_seeded_repair_lifecycle_and_missing_usage_are_honest(self):
        original = (self.root / "app.py").read_text()
        self.prepare()
        self.assertEqual(original, (self.root / "app.py").read_text())
        self.assertEqual("VALUE = 2\n", (self.workspace() / "app.py").read_text())
        self.assertTrue(tools.verify("trial")["seed_verified"])
        text = tools.read_context("trial", "app.py")
        (self.workspace() / "app.py").write_text("VALUE = 1\n")
        report = tools.grade("trial")
        self.assertTrue(report["passed"])
        self.assertEqual(len(text.encode()), report["instrumented_context_bytes"])
        self.assertIsNone(report["usage"]["input_tokens"])
        self.assertGreaterEqual(report["elapsed_seconds"], 0)
        self.assertEqual(original, (self.root / "app.py").read_text())

    def test_inherited_git_environment_cannot_modify_live_head_or_index(self):
        head = tools.git(self.root, "rev-parse", "HEAD")
        index = self.root / ".git/index"
        before = hashlib.sha256(index.read_bytes()).hexdigest()
        with patch.dict(os.environ, {"GIT_DIR": str(self.root / ".git"), "GIT_WORK_TREE": str(self.root), "GIT_INDEX_FILE": str(index)}):
            self.prepare()
            self.assertTrue(tools.verify("trial")["seed_verified"])
        self.assertEqual(head, tools.git(self.root, "rev-parse", "HEAD"))
        self.assertEqual(before, hashlib.sha256(index.read_bytes()).hexdigest())
        self.assertEqual("VALUE = 1\n", (self.root / "app.py").read_text())

    def test_rust_task_requires_the_specific_failing_test_header(self):
        with patch.object(tools, "SUITE", ROOT / ".agents/evals/tasks.json"):
            task = tools.tasks()["query-length-boundary"]
        misleading = "test hmm_alphabet_accepts_ambiguity_but_number_remains_strict ... ok\nother_test FAILED"
        self.assertFalse(all(marker in misleading for marker in task["expected_failure_contains"]))
        genuine = "---- hmm_alphabet_accepts_ambiguity_but_number_remains_strict stdout ----\nassertion failed"
        self.assertTrue(all(marker in genuine for marker in task["expected_failure_contains"]))

    def test_scope_violation_cannot_pass_even_when_tests_do(self):
        self.prepare()
        tools.verify("trial")
        (self.workspace() / "app.py").write_text("VALUE = 1\n")
        (self.workspace() / "extra.py").write_text("# unauthorized edit\n")
        report = tools.grade("trial")
        self.assertFalse(report["passed"])
        self.assertEqual(["extra.py"], report["scope_violations"])
        self.assertEqual(0, report["checks"][0]["returncode"])

    def test_agent_commits_do_not_hide_changed_paths(self):
        self.prepare()
        tools.verify("trial")
        (self.workspace() / "app.py").write_text("VALUE = 1\n")
        tools.git(self.workspace(), "add", ".")
        tools.git(self.workspace(), "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "--quiet", "-m", "repair")
        report = tools.grade("trial")
        self.assertTrue(report["passed"])
        self.assertEqual(["app.py"], report["changed_paths"])

    def test_unverified_or_unrelated_failures_do_not_establish_a_task(self):
        self.prepare()
        with self.assertRaisesRegex(ValueError, "verified"):
            tools.grade("trial")
        directory, state = tools.load_run("trial")
        state["task"]["expected_failure_contains"] = ["unrelated marker"]
        state["task_sha256"] = tools.digest(state["task"])
        tools.save_run(directory, state)
        self.assertFalse(tools.verify("trial")["seed_verified"])

    def test_changed_seed_is_rejected_before_verification(self):
        self.prepare()
        (self.workspace() / "app.py").write_text("VALUE = 1\n")
        with self.assertRaisesRegex(ValueError, "before editing"):
            tools.verify("trial")

    def test_usage_requires_actual_nonnegative_counts_and_provenance(self):
        file = self.root / "usage.json"
        for usage in ({"input_tokens": True, "output_tokens": 0, "source": "host"},
                      {"input_tokens": -1, "output_tokens": 2, "source": "host"},
                      {"input_tokens": 1, "output_tokens": 2, "source": ""}):
            file.write_text(json.dumps(usage))
            with self.assertRaises(ValueError):
                tools.usage_record(file)
        file.write_text(json.dumps({"input_tokens": 500, "output_tokens": 100, "source": "exported host usage"}))
        self.assertEqual(500, tools.usage_record(file)["input_tokens"])

    def test_paths_and_existing_runs_are_protected(self):
        for name in ("../outside", "UPPER", "", "x/y"):
            with self.assertRaises(ValueError):
                tools.run_directory(name)
        self.prepare()
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.prepare()
        with self.assertRaisesRegex(ValueError, "outside"):
            tools.local_path(self.workspace(), "../../outside")
        (self.workspace() / "link").symlink_to(self.root / "app.py")
        with self.assertRaisesRegex(ValueError, "outside"):
            tools.local_path(self.workspace(), "link")

    def test_prepare_does_not_export_dirty_live_work(self):
        (self.root / "app.py").write_text("VALUE = 999\n")
        self.prepare()
        self.assertEqual("VALUE = 2\n", (self.workspace() / "app.py").read_text())
        self.assertEqual("VALUE = 999\n", (self.root / "app.py").read_text())

    def test_seed_matching_fails_without_touching_live_source(self):
        task = copy.deepcopy(self.task)
        task["mutations"][0]["old"] = "not in source"
        self.suite.write_text(json.dumps({"schema_version": 1, "tasks": [task]}))
        with self.assertRaisesRegex(ValueError, "does not match"):
            self.prepare()
        self.assertEqual("VALUE = 1\n", (self.root / "app.py").read_text())

    def test_compare_rejects_unmatched_tasks_and_repeated_attempts(self):
        self.prepare("one")
        self.prepare("two", "navigation")
        reports = []
        for run in ("one", "two"):
            tools.verify(run)
            (self.workspace(run) / "app.py").write_text("VALUE = 1\n")
            tools.grade(run)
            reports.append(tools.run_directory(run) / "report-001.json")
        summary = tools.compare(reports)
        self.assertEqual({"baseline", "navigation"}, set(summary["variants"]))
        self.assertIsNone(summary["variants"]["baseline"]["median_input_tokens"])
        with self.assertRaisesRegex(ValueError, "one final attempt"):
            tools.compare([reports[0], reports[0]])
        value = json.loads(reports[1].read_text())
        value["baseline_commit"] = "different"
        reports[1].write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "same baseline"):
            tools.compare(reports)

    def test_timeout_is_recorded_and_not_a_valid_seed(self):
        self.prepare()
        directory, state = tools.load_run("trial")
        state["task"]["timeout_seconds"] = 1
        state["task"]["acceptance"] = [["{python}", "-c", "import time; time.sleep(10)"]]
        state["task_sha256"] = tools.digest(state["task"])
        tools.save_run(directory, state)
        result = tools.verify("trial")
        self.assertFalse(result["seed_verified"])
        self.assertTrue(result["checks"][0]["timed_out"])

    def test_suite_has_rust_and_tooling_tasks_without_paid_agents(self):
        with patch.object(tools, "SUITE", ROOT / ".agents/evals/tasks.json"):
            suite = tools.tasks()
        self.assertEqual(3, len(suite))
        self.assertEqual({"offline", "rust-offline"}, {task["tier"] for task in suite.values()})


if __name__ == "__main__":
    unittest.main()
