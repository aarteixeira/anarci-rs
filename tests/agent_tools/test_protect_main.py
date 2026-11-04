"""Branch-protection policy is explicit, authenticated and read-only by default."""
import contextlib
import copy
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("protect_main", ROOT / "scripts/protect_main.py")
tools = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tools)
POLICY = json.loads((ROOT / ".github/main-protection.json").read_text())


class ProtectionTests(unittest.TestCase):
    def test_policy_requires_pr_and_checks_even_for_admins(self):
        policy = POLICY["protection"]
        self.assertTrue(policy["enforce_admins"])
        self.assertTrue(policy["required_status_checks"]["strict"])
        self.assertIsNotNone(policy["required_pull_request_reviews"])
        self.assertFalse(policy["allow_force_pushes"])
        self.assertFalse(policy["allow_deletions"])
        self.assertEqual({"repository-map", "Fast validation"}, {c["context"] for c in policy["required_status_checks"]["checks"]})
        self.assertTrue(all(c["app_id"] == 15368 for c in policy["required_status_checks"]["checks"]))

    def test_api_shape_normalization_and_drift(self):
        desired = copy.deepcopy(POLICY["protection"])
        actual = copy.deepcopy(desired)
        for key in ("enforce_admins", "required_linear_history", "allow_force_pushes", "allow_deletions", "required_conversation_resolution"):
            actual[key] = {"enabled": actual[key]}
        actual["required_status_checks"]["checks"].reverse()
        self.assertEqual([], tools.differences(actual, desired))
        actual["enforce_admins"] = {"enabled": False}
        self.assertEqual(["enforce_admins"], tools.differences(actual, desired))
        self.assertEqual(["branch is not protected"], tools.differences(None, desired))

    def test_review_bypass_allowances_are_detected_as_drift(self):
        actual = copy.deepcopy(POLICY["protection"])
        actual["required_pull_request_reviews"]["bypass_pull_request_allowances"] = {"users": [{"login": "bypass-user"}], "teams": [], "apps": []}
        self.assertEqual(["required_pull_request_reviews"], tools.differences(actual, POLICY["protection"]))

    def test_default_check_never_mutates(self):
        with patch.object(tools, "api", return_value=POLICY["protection"]) as api, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(0, tools.main([]))
        api.assert_called_once_with("repos/aarteixeira/anarci-rs/branches/main/protection")

    def test_apply_writes_then_verifies(self):
        with patch.object(tools, "api", return_value=POLICY["protection"]) as api, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(0, tools.main(["--apply"]))
        self.assertEqual(2, api.call_count)
        self.assertEqual("PUT", api.call_args_list[0].args[1])
        self.assertEqual(1, len(api.call_args_list[1].args))

    def test_http_errors_are_not_silently_missing_protection(self):
        failure = subprocess.CompletedProcess([], 1, '{"status":"403","message":"Forbidden"}', 'gh: Forbidden')
        with patch.object(tools.subprocess, "run", return_value=failure), self.assertRaisesRegex(RuntimeError, "Forbidden"):
            tools.api("repos/test/repo/branches/main/protection")
        missing = subprocess.CompletedProcess([], 1, '{"status":"404","message":"Branch not protected"}', '')
        with patch.object(tools.subprocess, "run", return_value=missing):
            self.assertIsNone(tools.api("repos/test/repo/branches/main/protection"))


if __name__ == "__main__":
    unittest.main()
