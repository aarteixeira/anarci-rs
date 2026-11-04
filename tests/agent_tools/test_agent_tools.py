"""Offline contract tests for agent maps; no native build or installed package needed."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("agent_tools", ROOT / "scripts/agent_tools.py")
tools = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tools)


class RepositoryTests(unittest.TestCase):
    def test_no_tracked_files_match_ignore_rules(self):
        result = subprocess.run(["git", "ls-files", "--cached", "--ignored", "--exclude-standard", "-z"], cwd=ROOT, check=True, stdout=subprocess.PIPE)
        self.assertEqual(b"", result.stdout, "Untrack local artifacts with git rm --cached; preserve their local copies")

    def test_distribution_ignore_boundaries(self):
        excluded = (
            "workbench/pdb-ab-census/analysis.py", "outputs/release/results.json",
            "output/figure.pdf", "reports/report.md", "target/build.rs", "dist/package.whl",
            "PLAN.md", "LEARNED_ENGINE_NOTEBOOK.md", "docs/PYPI_RELEASE_PLAN.md",
            "docs/RELEASE_EVIDENCE.md", "docs/PERFORMANCE.md", "scripts/anarcii_baseline.py",
            "scripts/measure_partial_recall.py", "scripts/train_gnode2.sbatch",
            "reference_data/ground_truth/imgt.json", "reference_data/dat/HMMs/ALL.hmm.h3m",
            ".DS_Store", ".env", ".venv/bin/python", "checkpoints/model.safetensors",
        )
        retained = (
            "Cargo.lock", "Cargo.toml", "pyproject.toml", "AGENTS.md", ".agents/repo.json",
            "docs/agents/repo-map.json", "docs/DISTRIBUTION.md", "docs/FEATURE_COVERAGE.md",
            "reference_data/dat/HMMs/ALL.hmm", "reference_data/dat/HMMs/FEW.hmm",
            "reference_data/germlines.py", "crates/anarci-core/data/germlines.json",
            "tests/fixtures/golden.json.gz", "tests/python/safe_expected.json",
            "scripts/benchmark_release.py", "scripts/gen_constants.py", "licenses/HMMER-LICENSE.txt",
        )
        for paths, expected in ((excluded, 0), (retained, 1)):
            for path in paths:
                with self.subTest(path=path):
                    result = subprocess.run(["git", "check-ignore", "--no-index", "--quiet", "--", path], cwd=ROOT)
                    self.assertEqual(expected, result.returncode, path)

    def test_graph_references_are_valid(self):
        graph = tools.build_graph(ROOT)
        ids = {node["id"] for node in graph["nodes"]}
        self.assertEqual(len(ids), len(graph["nodes"]))
        for edge in graph["edges"]:
            self.assertIn(edge["from"], ids)
            self.assertIn(edge["to"], ids)
            self.assertTrue((ROOT / edge["evidence"]).is_file(), edge)
        self.assertTrue(any(edge["kind"] == "imports" for edge in graph["edges"]))
        self.assertTrue(any(edge["kind"] == "flow" for edge in graph["edges"]))

    def test_deterministic_rendering(self):
        self.assertEqual(tools.rendered(ROOT), tools.rendered(ROOT))

    def test_core_change_reaches_all_consumers(self):
        result = tools.impact(ROOT, ["crates/anarci-core/src/orchestrate.rs"])
        self.assertTrue({"anarci-core", "anarci-hmm", "anarci-cli", "anarci-py", "anarci-learned", "python-api"} <= set(result["affected_components"]))
        self.assertIn("core", result["suggested_recipes"])
        self.assertIn("python", result["suggested_recipes"])
        self.assertEqual([], result["unmapped_paths"])

    def test_profile_changes_reach_embedded_consumers(self):
        for profile in ("ALL", "FEW"):
            result = tools.impact(ROOT, [f"reference_data/dat/HMMs/{profile}.hmm"])
            self.assertTrue({"reference-data", "anarci-cli", "anarci-hmm", "anarci-py", "python-api"} <= set(result["affected_components"]))
            self.assertIn("cli", result["suggested_recipes"])

    def test_context_is_bounded_and_reports_neighbors(self):
        result = tools.context(ROOT, "anarci-core", limit=2)
        self.assertEqual(2, len(result["files"]))
        self.assertTrue(result["truncated"])
        self.assertTrue(any(edge["from"] == "anarci-hmm" for edge in result["consumers"]))
        result = tools.context(ROOT, "python/anarci_rs/_api.py")
        self.assertTrue(result["file_imports"])

    def test_deleted_file_and_unmapped_path(self):
        result = tools.impact(ROOT, ["crates/anarci-cli/src/removed.rs", "unknown/example.txt"])
        self.assertIn("anarci-cli", result["affected_components"])
        self.assertEqual(["unknown/example.txt"], result["unmapped_paths"])

    def test_release_manifest_is_conservative(self):
        result = tools.impact(ROOT, ["Cargo.lock"])
        self.assertTrue({"release", "hmmer-sys", "anarci-learned", "python-api"} <= set(result["affected_components"]))

    def test_agent_routes_do_not_trigger_native_build_recipes(self):
        result = tools.impact(ROOT, ["scripts/agent_tools.py", "docs/agents/PLAYBOOK.md", "AGENTS.md", ".github/workflows/agents.yml"])
        self.assertEqual(["agent-infra"], result["affected_components"])
        self.assertEqual(["agent"], list(result["suggested_recipes"]))

    def test_outside_paths_rejected(self):
        with self.assertRaisesRegex(ValueError, "outside repository"):
            tools.impact(ROOT, ["../outside.rs"])

    def test_impact_never_runs_commands(self):
        with patch.object(tools, "git_inventory", return_value=None), patch.object(tools.subprocess, "run", side_effect=AssertionError("validation command executed")):
            tools.impact(ROOT, ["crates/hmmer-sys/build.rs"])

    def test_inventory_omits_binary_fixtures_and_artifacts(self):
        files = set(tools.source_files(ROOT))
        self.assertNotIn("tests/fixtures/golden.json.gz", files)
        self.assertFalse(any(path.startswith("target/") or "__pycache__" in path for path in files))
        self.assertFalse(set(tools.OUTPUTS) & files)


class SyntheticRepositoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.put("Cargo.toml", '[workspace]\nmembers = ["crates/a", "crates/b"]\n')
        self.put("crates/a/Cargo.toml", '[package]\nname = "a"\n[dependencies]\nb = { path = "../b", optional = true }\n[features]\nextra = ["dep:b"]\n')
        self.put("crates/b/Cargo.toml", '[package]\nname = "b"\n')
        self.put("python/pkg/__init__.py", "from . import api\n")
        self.put("python/pkg/api.py", "from . import helper\n")
        self.put("python/pkg/helper.py", "VALUE = 1\n")
        self.put(".agents/repo.json", json.dumps({"schema_version": 1, "components": [{"id": "python", "path": "python", "description": "Python", "depends_on": ["a"], "recipes": []}], "crate_recipes": {"a": ["a-test"], "b": []}, "recipes": {"a-test": {"commands": ["cargo test -p a"], "requires": "Rust", "tier": "targeted"}}, "routes": [], "flows": []}))

    def put(self, path, text):
        file = self.root / path
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text(text)

    def test_optional_dependency_and_relative_imports(self):
        graph = tools.build_graph(self.root)
        self.assertTrue(any(e["from"] == "a" and e["to"] == "b" and e.get("optional") for e in graph["edges"]))
        self.assertTrue(any(e["kind"] == "imports" and e["from"] == "file:python/pkg/api.py" and e["to"] == "file:python/pkg/helper.py" for e in graph["edges"]))
        result = tools.impact(self.root, ["crates/b/src/deleted.rs"])
        self.assertEqual(["a", "b", "python"], result["affected_components"])

    def test_check_detects_new_file_and_recovers(self):
        with patch.object(tools, "ROOT", self.root), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(1, tools.main(["check"]))
            self.assertEqual(0, tools.main(["generate"]))
            self.assertEqual(0, tools.main(["check"]))
            self.put("crates/a/src/new.rs", "pub fn new() {}\n")
            self.assertEqual(1, tools.main(["check"]))
            self.assertEqual(0, tools.main(["generate"]))
            self.assertEqual(0, tools.main(["check"]))

    def test_source_body_changes_do_not_invalidate_topology(self):
        before = tools.rendered(self.root)
        self.put("python/pkg/helper.py", "VALUE = 2\n")
        self.assertEqual(before, tools.rendered(self.root))
        self.put("python/pkg/helper.py", "from . import api\n")
        self.assertNotEqual(before, tools.rendered(self.root))

    def test_invalid_route_fails(self):
        config = json.loads((self.root / ".agents/repo.json").read_text())
        config["routes"] = [{"patterns": ["*"], "components": ["missing"]}]
        self.put(".agents/repo.json", json.dumps(config))
        with self.assertRaisesRegex(ValueError, "Unknown route component"):
            tools.build_graph(self.root)

    def test_unknown_local_dependency_fails(self):
        self.put("crates/a/Cargo.toml", '[package]\nname = "a"\n[dependencies]\nmissing = {path = "../missing"}\n')
        with self.assertRaisesRegex(ValueError, "Unmapped local dependency"):
            tools.build_graph(self.root)

    def test_git_ignored_source_does_not_change_maps(self):
        subprocess.run(["git", "init"], cwd=self.root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.put(".gitignore", "scripts/local-cache.json\npython/pkg/ignored.py\n")
        before = tools.rendered(self.root)
        self.put("scripts/local-cache.json", '{"private": true}')
        self.put("python/pkg/ignored.py", "not even valid python!")
        self.put("local-notes.md", "Private scratch notes")
        self.assertEqual(before, tools.rendered(self.root))

    def test_symlinked_base_directory_is_not_read(self):
        with tempfile.TemporaryDirectory() as outside:
            (Path(outside) / "private.py").write_text("broken python!")
            (self.root / "scripts").symlink_to(outside, target_is_directory=True)
            self.assertNotIn("scripts/private.py", set(tools.source_files(self.root)))

    def test_symlinked_source_is_not_read(self):
        (self.root / "python/pkg/outside.py").symlink_to("/does/not/exist")
        self.assertNotIn("python/pkg/outside.py", set(tools.source_files(self.root)))

    def test_git_collects_staged_unstaged_and_untracked_paths_with_spaces(self):
        def git(*args):
            subprocess.run(["git", "-c", "core.hooksPath=/dev/null", "-c", "commit.gpgsign=false", *args], cwd=self.root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        git("init")
        git("config", "user.email", "agent-test@example.invalid")
        git("config", "user.name", "Agent test")
        self.put("tracked.txt", "initial")
        self.put("deleted.txt", "initial")
        self.put(".gitignore", "ignored.txt\n")
        git("add", ".")
        git("commit", "-m", "fixture")
        self.put("tracked.txt", "changed")
        self.put("staged file.txt", "staged")
        git("add", "staged file.txt")
        self.assertEqual(self.root / "deleted.txt", (self.root / "deleted.txt").resolve())
        git("rm", "deleted.txt")
        self.put("untracked file.txt", "untracked")
        self.put("ignored.txt", "ignored")
        (self.root / "scripts").mkdir()
        self.assertEqual(self.root / "python/pkg/helper.py", (self.root / "python/pkg/helper.py").resolve())
        self.assertEqual(self.root / "scripts/moved.py", (self.root / "scripts/moved.py").resolve())
        git("mv", "python/pkg/helper.py", "scripts/moved.py")
        self.assertEqual(["deleted.txt", "python/pkg/helper.py", "scripts/moved.py", "staged file.txt", "tracked.txt", "untracked file.txt"], tools.git_paths(self.root))


if __name__ == "__main__":
    unittest.main()
