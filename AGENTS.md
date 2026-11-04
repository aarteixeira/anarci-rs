# Agent development guide

Make the smallest change that satisfies the task. Read the affected sources and record the existing worktree changes before editing; preserve work you did not create. Do not reset, clean, reformat unrelated code, or commit without instruction.

## Start here

- Use the [playbook](docs/agents/PLAYBOOK.md) for routing and validation, and the [task template](docs/agents/TASK_TEMPLATE.md) for scoped work and handoffs.
- Read the generated [architecture](docs/agents/ARCHITECTURE.md), then use `context` to retrieve the files and dependencies for the component you need. The full [repository map](docs/agents/repo-map.json) is for programmatic queries; avoid loading it all into context. Verify details in source. [Repository configuration](.agents/repo.json) holds curated recipes and boundaries.
- Run the stdlib-only [agent helper](scripts/agent_tools.py) from the repository root with Python 3.11+. `python3` is an alternative only if its version is at least 3.11. This tooling requirement does not change the package's supported Python versions.

```sh
python3.11 scripts/agent_tools.py doctor
python3.11 scripts/agent_tools.py context anarci-core
python3.11 scripts/agent_tools.py impact crates/anarci-core/src/orchestrate.rs
python3.11 scripts/agent_tools.py impact --git
python3.11 scripts/agent_tools.py generate
python3.11 scripts/agent_tools.py check
```

`generate` refreshes the two deterministic generated documents; edit their inputs, not their output. `check` detects drift. `impact [paths...]` or `impact --git` suggests tests only: it never runs them. Its conservative crate reverse dependencies are not a full call graph. Inspect suggestions before choosing validation.

## Invariants

- Keep core numbering independent of native FFI. Preserve the in-process HMMER backend, verified source acquisition, native ownership and thread-safety rules, and input validation before C calls. Native failure must remain distinct from a successful no-hit result.
- Preserve the separate `anarci` compatibility and `anarci_rs` extended contracts, including defaults, tuple shapes, ordering, errors and output formats. Compare exact/ALL mode with the pinned reference; pan mode is a separate option. See [compatibility](docs/COMPATIBILITY.md).
- Use identified public or synthetic validation inputs and record provenance, settings and tolerances. Do not replace expected results merely to pass tests. A passing fixture suite does not prove general biological accuracy; geometry does not establish binding specificity.
- Keep learned/model features optional and experimental. No validated trained checkpoint ships with the package. Synthetic codec tests or training accuracy do not establish held-out numbering accuracy. See [experimental tools](docs/EXPERIMENTAL.md).

## Permission and scope

Do not publish, push, dispatch remote CI, install dependencies, download data, start training/GPU jobs, or run costly benchmarks without authorization for that work. Native builds may download HMMER even with Cargo offline: inspect [build prerequisites](docs/BUILDING.md) and use a verified local `HMMER_TARBALL` and cached dependencies when network access is not authorized. Do not run `cargo clean` as routine troubleshooting.

Keep local research workbenches, analysis reports, run outputs and historical planning notes out of version control under the [repository distribution policy](docs/DISTRIBUTION.md). Retain runtime resources, licenses and reusable regression fixtures; they are source dependencies, not disposable outputs. Treat historical fixtures as separate from release evidence. Do not inspect or package unrelated data. Preserve attribution and artifact exclusions; publication requires approval for the exact artifacts and index under the [release procedure](docs/RELEASING.md).

For parallel work, assign disjoint file ownership and explicit interface assumptions. Each worker reports changes, checks and unresolved risks; one integrator reviews shared boundaries and runs the agreed integration checks. Never overwrite another worker's edits to resolve a conflict.
