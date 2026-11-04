# Task brief and handoff

Use with the [root guide](../../AGENTS.md), [playbook](PLAYBOOK.md), [architecture](ARCHITECTURE.md) and [curated recipes](../../.agents/repo.json). Keep only fields relevant to the task; fill evidence with observed results, not planned success.

## Brief

- Objective and user-visible behavior:
- Acceptance conditions (specific inputs, expected outputs or failure behavior):
- Owned files and shared interfaces; other workers and integrator:
- Existing worktree changes to preserve:
- Sources read and current behavior confirmed:
- In scope:
- Out of scope and permission boundaries (network, installs, publication, remote CI, datasets, training, benchmark budget):

## Contracts and evidence plan

- Native invariants affected: validation before C, ownership/thread safety, failure versus no hit, verified HMMER source:
- Python/API invariants affected: compatibility versus extended defaults, shapes, ordering, exceptions, supported interpreters:
- Scientific claim and limits: input provenance, pinned reference/profile hashes, exact versus pan settings, comparison tolerances:
- Experimental boundary: optional features, split leakage checks if relevant, no unsupported trained-model or biological accuracy claim:
- Impact command and suggestions reviewed (suggestions do not run tests; crate reverse dependencies are conservative, not a full call graph):
- Validation tiers selected and why; prerequisites and time/resource limit:
- Acceptance evidence required before completion:

Run helper commands from the root using Python 3.11+; `python3` is acceptable only after verifying its version.

```sh
python3.11 scripts/agent_tools.py impact path/to/changed_file
python3.11 scripts/agent_tools.py impact --git
python3.11 scripts/agent_tools.py check
```

## Completion and handoff

- Files created/changed and behavior implemented:
- Interface assumptions or coordination needed by the next worker:
- Checks actually run: exact command, environment/features, exit status and result:
- Evidence locations (test output, fixture provenance, artifact hashes when applicable):
- Checks skipped/blocked and reason; no implied pass:
- Acceptance conditions met or still open:
- Remaining risks, scientific limits and next action/owner:
- Combined-diff review and integration checks still required:
