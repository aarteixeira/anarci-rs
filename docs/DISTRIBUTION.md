# Repository distribution policy

Version the software and the resources needed to build, test and maintain it.
Keep local research and generated results out of the repository. The
[ignore rules](../.gitignore) enforce this boundary; Python artifacts have a
separate allowlist in the [build manifest](../pyproject.toml) and
[release validator](../scripts/validate_release.py).

## Keep in version control

- Rust and Python source, portable build configuration, CI and the lockfile.
- Runtime HMM source profiles, germline inputs, generated runtime constants and
  their licenses and attribution. Generated files that the build consumes remain
  source dependencies until the build can reproduce them independently.
- Reusable tests, their regression fixtures, generators and validation tools.
  Historical scientific fixtures remain available for opt-in development tests;
  they are excluded from Python release artifacts and are not release evidence.
- API/build documentation, agent instructions, repository maps and task templates.
  The maps are intentional development artifacts, regenerated and checked in CI.

## Keep local

Research workbenches, PDB/patent/reagent analyses, reports and figures; release
run logs, benchmark measurements and artifact manifests; historical plans and
notebooks; paper-specific experiments and workstation launch notes; build
products, caches, environments and model checkpoints. Avoid saving new analysis
outputs beside source files: use the ignored output directories instead.

The old ground-truth capture directory is local analysis data. The regression
fixtures used by the Rust suites remain versioned separately. Do not apply broad
ignore rules to JSON, Markdown, compressed data or all reference data: those
formats also contain required source and test inputs.

## HMMER indexes

The native engine embeds the plain ALL and FEW HMM profiles. Pressed HMMER
indexes are regenerable and ignored. External-reference tools may need them;
for example, before running the pan comparison tool in a suitable isolated
reference environment:

```sh
hmmpress reference_data/dat/HMMs/FEW.hmm
python scripts/validate_pan.py
```

Neither command is part of routine agent-map validation. The first requires an
installed external HMMER toolchain; the second uses historical regression inputs.
Review prerequisites and input provenance before running either.

## Before committing

Review staged paths and run the offline boundary tests:

```sh
git diff --cached --stat
git ls-files --cached --ignored --exclude-standard
python3.12 -m unittest discover -s tests/agent_tools -v
python3.12 -m unittest discover -s tests/packaging -v
python3.12 scripts/agent_tools.py check
```

Use Python 3.11+ for agent tooling. The tracked/ignored query should print nothing.
For a file that is already tracked, ignore rules alone do not stop distribution:
remove it from the index with `git rm --cached` after verifying the path, preserving
the local copy. Regenerate repository maps after changing the boundary.

Untracking files changes future checkouts, not earlier commits. Do not rewrite
published history without a separate review and explicit authorization. If
sensitive material was ever committed, removing it from the current tree does
not revoke access to it; handle that as a separate incident.
