#!/usr/bin/env python3
"""Check managed GitHub branch protections; mutation requires explicit --apply."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / ".github/main-protection.json"


def api(endpoint, method="GET", payload=None):
    command = ["gh", "api", "--method", method, endpoint]
    if payload is not None:
        command.extend(["--input", "-"])
    result = subprocess.run(command, input=None if payload is None else json.dumps(payload), text=True, capture_output=True)
    if result.returncode:
        try:
            error = json.loads(result.stdout)
        except ValueError:
            error = {}
        if method == "GET" and str(error.get("status")) == "404" and error.get("message") == "Branch not protected":
            return None
        detail = result.stderr.strip() or error.get("message", "gh api failed")
        if error.get("errors"):
            detail += ": " + json.dumps(error["errors"])
        raise RuntimeError(detail)
    return json.loads(result.stdout)


def normalized(actual):
    if actual is None:
        return None
    checks = actual.get("required_status_checks")
    reviews = actual.get("required_pull_request_reviews")
    keys = ("dismiss_stale_reviews", "require_code_owner_reviews", "required_approving_review_count", "require_last_push_approval")
    normalized_reviews = None if reviews is None else {key: reviews.get(key, False) for key in keys}
    if normalized_reviews is not None:
        allowances = reviews.get("bypass_pull_request_allowances", {})
        normalized_reviews["bypass_pull_request_allowances"] = {
            kind: sorted(str(actor.get("login", actor.get("slug", actor.get("id")))) if isinstance(actor, dict) else str(actor) for actor in allowances.get(kind, []))
            for kind in ("users", "teams", "apps")
        }
    result = {
        "required_status_checks": None if checks is None else {
            "strict": checks["strict"],
            "checks": sorted(checks.get("checks", []), key=lambda item: (item["context"], item.get("app_id", -1))),
        },
        "required_pull_request_reviews": normalized_reviews,
        "restrictions": actual.get("restrictions"),
    }
    for key in ("enforce_admins", "required_linear_history", "allow_force_pushes", "allow_deletions", "required_conversation_resolution"):
        value = actual.get(key, False)
        result[key] = value.get("enabled", False) if isinstance(value, dict) else value
    return result


def differences(actual, desired):
    observed = normalized(actual)
    expected = normalized(desired)
    if observed is None:
        return ["branch is not protected"]
    return [key for key in expected if observed[key] != expected[key]]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Read-only drift check (default)")
    mode.add_argument("--apply", action="store_true", help="Apply the reviewed policy to the named repository and branch")
    args = parser.parse_args(argv)
    try:
        policy = json.loads(POLICY.read_text())
        endpoint = f"repos/{policy['repository']}/branches/{policy['branch']}/protection"
        if args.apply:
            api(endpoint, "PUT", policy["protection"])
        drift = differences(api(endpoint), policy["protection"])
        print(json.dumps({"repository": policy["repository"], "branch": policy["branch"], "applied": args.apply, "drift": drift}, indent=2))
        return 1 if drift else 0
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(f"protect_main: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
