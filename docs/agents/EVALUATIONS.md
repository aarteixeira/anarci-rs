# Evaluate agent development tasks

Use the [evaluation runner](../../scripts/agent_eval.py) to compare repair correctness,
retrieved context and wall time on the same source revision. It creates disposable
source snapshots under the ignored output directory. It does not launch an agent,
call a model API, change the live checkout or publish results.

## Run one trial

Commit the intended task definitions and source first: snapshots come from a pinned
Git commit, not dirty local edits. List the tasks, then prepare and verify a seed:

```sh
python3.12 scripts/agent_eval.py list
python3.12 scripts/agent_eval.py prepare impact-consumers trial-navigation-01 \
  --variant navigation --ref HEAD --agent your-agent-version --model your-model-id
python3.12 scripts/agent_eval.py verify trial-navigation-01
```

`prepare` exports committed source, introduces one controlled defect and creates a
local Git baseline containing that defect. It refuses to overwrite an existing run.
The [task suite](../../.agents/evals/tasks.json) defines the brief, allowed edit paths,
acceptance commands and expected failure. If the seed no longer matches the source,
preparation fails rather than silently benchmarking a different task.

`verify` must observe the intended failing test before timing starts. A missing
compiler, dependency or unrelated failure is not a valid seed. The Rust task uses
only the offline core contract target; the other tasks need Python 3.11+ and Git.
No biological fixture captures, native HMMER build or model training are required.

Give the agent the emitted workspace path and task brief. Keep it inside that
workspace, with only the allowed source edits. For measured retrieval, use:

```sh
python3.12 scripts/agent_eval.py read trial-navigation-01 scripts/agent_tools.py \
  --offset 210 --limit 70
python3.12 scripts/agent_eval.py context trial-navigation-01 anarci-core --limit 10
```

The `baseline` arm should use ordinary source reads; the `navigation` arm can use
the repository maps and symbol helpers. Use the same retrieval policy across
repeats. Reads performed outside the evaluator are not counted in its byte metric.
Run the evaluator from the original checkout, not the agent's snapshot.

After the repair:

```sh
python3.12 scripts/agent_eval.py grade trial-navigation-01
```

Grading runs bounded acceptance commands with per-command logs and checks changed
paths against the task scope, including changes committed by the agent. A passing
run must make a source change, pass every acceptance check and stay in scope.
Python cache paths are isolated per check so a rapid same-size source edit cannot
reuse the seeded bytecode. Inherited Git repository-selection variables are removed
from snapshot and acceptance commands to protect the original index and branch.

## Record usage and compare trials

The runner cannot observe the model's token usage. Missing counts stay null; source
bytes are never converted into pretend token counts. If the agent host exports
usage, supply a JSON object with nonnegative integer `input_tokens`, `output_tokens`
and a `source` string identifying that export:

```sh
# Use this instead of the earlier grade command when usage is available.
python3.12 scripts/agent_eval.py grade trial-navigation-01 \
  --usage outputs/agent-evals/exported-usage.json
python3.12 scripts/agent_eval.py compare \
  outputs/agent-evals/trial-baseline-01/report-001.json \
  outputs/agent-evals/trial-navigation-01/report-001.json
```

Use one final attempt per independent run. Comparisons require the same baseline
commit, task definition and agent/model/settings metadata. Pass `--settings` to
`prepare` with a JSON object recording controlled settings such as seed, reasoning
effort and token budget. Unspecified metadata is visibly unspecified, not inferred.

Reports include correctness, scope violations, wall time since seed verification,
bytes returned by instrumented retrieval, optional exported token counts and the
exact check commands. Wall time includes idle time and ends when grading starts;
it is not CPU time or active thinking time. The comparator reports counts and
medians without claiming causality or statistical significance.

Repeat paired trials across all tasks with the same agent version, model, budget,
prerequisites and baseline revision. Alternate arm order. Report failures as well
as successes, and inspect patches for test-specific shortcuts. Three controlled
repair tasks are a starting benchmark, not a measure of autonomous research or
production engineering competence. The task seeds and tests are visible; this is
not a blind benchmark or an adversarial security sandbox.

Keep raw workspaces, transcripts, usage exports and reports local. Version the
reusable suite and evaluator, not run results, under the
[distribution policy](../DISTRIBUTION.md).
