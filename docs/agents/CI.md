# CI gates and protected main

Every PR must pass `repository-map` and `Fast validation` before merging into
`main`. Both checks are bound to the GitHub Actions app, and the branch must be up
to date. The [managed protection policy](../../.github/main-protection.json) also
requires a PR, enforces the rules for admins, blocks force pushes and deletion,
requires resolved conversations and keeps history linear.

The approval count is zero because the repository has a single maintainer; requiring
that maintainer to approve their own PR would block normal work. This does not remove
the PR or CI requirements. User/team/app PR-bypass allowances must remain empty.
Review changes to tests, CI and the protection policy carefully; checks that a PR
can edit are not a security boundary against a malicious contributor.

## Fast PR checks

[Agent repository contracts](../../.github/workflows/agents.yml) run on every PR:
repository maps, distribution boundaries, source archive safety and offline tests
for navigation, evaluation and protection tooling. No path filter may skip the
required check.

[Fast source validation](../../.github/workflows/fast-validation.yml) runs:

- Existing core unit/input contracts and generated core/learned property contracts,
  without historical biological fixture captures or model features.
- Native HMMER robustness and shared CLI library contracts on Linux.
- One CPython 3.12 development wheel, installed and tested outside the checkout.

`Fast validation` is an always-run aggregate job. It succeeds only if both Rust
and native/Python jobs succeed; skipped, cancelled or failed dependencies do not
count as passing. Jobs have timeouts and cache Cargo/build outputs. A newer push
cancels superseded fast runs, not the current head's required checks.

This is a smaller required suite, not a promise that a fresh native build is
instant. Optional model and parser tests retain their documented feature/dependency
requirements; a default build does not qualify those optional paths.

## Full qualification

The [candidate distribution matrix](../../.github/workflows/wheels.yml) runs on
relevant pushes to `main`, or by explicit workflow dispatch. It covers all configured
wheel platforms/interpreters, independent source installation and optional CPU model
bindings. It no longer starts the entire matrix for every PR.

A green fast gate does not establish platform qualification or publication readiness.
Review the completed full matrix and exact artifacts under the
[release procedure](../RELEASING.md) before publishing. A cancelled platform job is
not a passing platform result. No workflow publishes to a package index.

## Check or apply GitHub settings

These commands need an authenticated `gh` client with repository access:

```sh
python3.12 scripts/protect_main.py --check
# Only after explicit authorization for the settings change:
python3.12 scripts/protect_main.py --apply
```

The [helper](../../scripts/protect_main.py) is read-only by default, reports drift
and does not weaken protections to make a merge succeed. `--apply` writes the
reviewed policy and reads it back before reporting success. Authentication and
permission errors remain errors; they are not treated as an unprotected branch.
The settings live on GitHub, so editing the JSON alone does not change them.
