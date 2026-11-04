# Development playbook

Follow the [root guide](../../AGENTS.md). Start with an acceptance condition, owned paths and a validation budget in the [task template](TASK_TEMPLATE.md). Inspect the [architecture](ARCHITECTURE.md), [map](repo-map.json) and [curated recipes](../../.agents/repo.json), then read the code you will change.

## Route the task

| Change | Start with | Check across the boundary |
|---|---|---|
| Schemes, states, germlines, orchestration | [anarci-core](../../crates/anarci-core/) | Exact numbering, coordinates, insertions, species filters, batch ordering |
| Native scans or build | [anarci-hmm](../../crates/anarci-hmm/), [hmmer-sys](../../crates/hmmer-sys/) | C ownership, per-thread mutable state, input rejection, error propagation, verified tarball |
| Python API or compatibility | [Python packages](../../python/), [anarci-py](../../crates/anarci-py/) | Defaults, return shapes, exceptions, namespace isolation, installed extension |
| Files, DNA, structures, output or CLI | [anarci-cli](../../crates/anarci-cli/), [Python CLI](../../python/anarci_rs/cli.py) | Shared Rust/Python behavior, coordinates, output schema, malformed input |
| Codec, decoder, datasets or model | [anarci-learned](../../crates/anarci-learned/), [experimental bindings](../../crates/anarci-py/src/experimental.rs) | State conventions, deterministic splits, optional features, no production accuracy claim |
| Packaging or release | [build guide](../BUILDING.md), [release procedure](../RELEASING.md), [manifest](../../pyproject.toml) | Artifact allowlists, licenses, clean installation, extracted source build |

For cross-cutting changes, trace callers and consumers yourself. Run the [helper](../../scripts/agent_tools.py) from the root:

```sh
python3.11 scripts/agent_tools.py impact path/to/changed_file
python3.11 scripts/agent_tools.py impact --git
```

Use Python 3.11+ (`python3` is acceptable after checking its version). Impact reports conservative crate reverse dependencies and suggested tests, never executes them, and does not resolve a full call graph. Unknown paths or incomplete coverage require manual review, not a conclusion that no tests apply.

Use `context anarci-core` (or a source path) for a bounded file list, dependencies and consumers. The full map is for programmatic queries. Freshness checks cover file inventory, declared dependencies, static Python imports and curated policy; they do not validate behavior or automatically detect changes to reviewed flow semantics. Update curated flows, data dependencies and recipes when those relationships change. Commit intended new source files together with regenerated maps; ignored local artifacts are excluded.

## Validation tiers

Choose checks by changed behavior and available prerequisites. Commands below are examples, not an instruction to run every suite. Check curated recipes and input provenance first; avoid blanket workspace tests that may use historical fixtures.

1. **Static/local:** review the diff, links and scope; run `python3.11 scripts/agent_tools.py check`. After changing map inputs, run `generate` and then `check`; review both generated outputs. Use `doctor` to identify missing prerequisites, not to install them.
2. **Focused contracts:** run the smallest relevant synthetic or approved unit tests. Examples: `cargo test --locked -p anarci-core --lib`, `cargo test --locked -p anarci-core --test native_input_contract`, `cargo test --locked -p anarci-hmm --lib robustness_tests`, and `cargo test --locked -p anarci-cli --lib`. Native commands can trigger downloads or compilation; verify caches, local HMMER source and budget first.
3. **Integration/scientific:** test an installed extension with `python3.11 -m unittest discover -s tests/python -v` in the intended environment. For compatibility changes, use the pinned [reference settings and provenance](../reference-version.json), exact/ALL mode and documented tolerances. Record failures and limits; do not generalize fixture agreement to untested sequences. Optional model checks require their own feature build and approval for any costly work.
4. **Artifact/release:** follow the release procedure for clean wheel installs and builds from extracted source, archive inspection, attribution, and platform evidence. Record actual artifact hashes and environment details with the reviewed release artifacts only after checks finish; keep local run reports out of version control. Local success does not qualify other platforms. Upload remains a separate authorization.

If a check is blocked, report the missing prerequisite and the command not run. Do not mark it passed. Benchmark only with an agreed dataset and budget; compare equal profile sets, settings, workers and duplicate rates.

## Parallel work and handoff

Split independent changes by owned files. Agree on signatures, data shapes and acceptance tests before editing a shared interface. Keep one owner for manifests, generated outputs and integration; coordinate rather than letting several workers regenerate or reformat the same files.

Each handoff includes owned paths, behavior changed, interface assumptions, exact checks and outcomes, skipped checks with reasons, and remaining risks. The integrator reviews the combined diff against the original worktree, resolves interface mismatches without discarding others' edits, and reruns checks affected by the merge.
