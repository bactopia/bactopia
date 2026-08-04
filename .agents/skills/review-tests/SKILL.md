---
name: review-tests
description: Review nf-test run results and present a diagnostic summary with grouped error analysis. Use when asked to review tests, check test results, show test failures, analyze test output, investigate why tests failed, see what's broken, or check test status. Runs are multi-profile (docker/conda/singularity) with docker as the reference baseline. Accepts an optional timestamp argument to review a specific run.
---

# Review Tests

Run the review-tests CLI and present the results to the user.

Runs are **multi-profile**: each component/tier is tested across up to four profiles
-- `docker` (the reference baseline), `conda`, `singularity_galaxy`, `singularity_pull`.
A cell is one (component, tier, profile). Most interpretation is a comparison of each
profile against docker.

## Steps

1. Run `bactopia-review-tests` via the wrapper script using the **default text output**
   (do NOT use `--json`):
   ```
   bash .agents/skills/review-tests/scripts/run-bactopia-review-tests.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --silent
   ```
   If the user provided a timestamp argument (e.g., `/review-tests 20260324_081306`),
   add `--run 20260324_081306`.

2. Present the text output directly to the user. The CLI produces a clean summary with a
   **"Status Breakdown by Profile"** matrix and one section per failing status. Do NOT
   parse JSON or write extra code to reformat -- just relay the output with your interpretation.

3. Add interpretation and context after showing the output. Interpret **by status**
   (see the status reference below), and always frame failures as "which profiles differ
   from docker, and why." Summarize actionable items and next steps.

4. If the text output is too large for a single response, summarize the key sections
   (overview, status-by-profile matrix, failures) and note that per-file detail is in
   `summary.json`. Use `--json`/`--pretty` or read `summary.json` directly for structured detail.

## Status reference (per-cell `status`)

- **`passed`** -- cell matched the committed snapshot / assertions.
- **`version_drift` / `output_drift` / `version+output_drift`** -- this profile's outputs differ
  from docker: `version_drift` = `versions.yml` (runtime resolved a different tool version than the
  docker container pin), `output_drift` = tool output file(s), `version+output_drift` = both.
  `reason` names the divergent fields. When **docker passed** and only conda/singularity drift, this
  is **genuine dependency-solve divergence**, not a bug. Fix = the **sccmec-pattern test migration**
  (md5 the profile-stable files, existence-check the divergent ones, `versions` -> `contains('<tool>')`)
  -- **NOT** snapshot regeneration. See `files[]` / `suggested_edit` in `summary.json` for the exact
  bucketing. (A pure `version_drift` may instead warrant updating the container version pin.)
- **`snapshot_mismatch`** -- the snapshot didn't match but the files matrix could not attribute it to
  specific fields (drift not subclassified). Inspect `files[]` and `{profile}/stderr.txt`.
- **`snapshot_stale`** -- the committed `.snap` no longer matches the **reference runtime (docker)**;
  it shows on **all profiles including docker**. Fix: re-run with `--generate` under docker to
  re-baseline. NOT a content or tool problem (on `generate=false`, docker's own mismatch is promoted
  to this). Typical cause: a test's `snapshot()` shape was edited but `.snap` was never regenerated.
- **`assertion_failed`** -- a non-snapshot assertion failed (no output divergence detected). A test
  logic/assertion issue, not drift -- read `{profile}/stderr.txt`.
- **`non_reproducible`** -- two docker runs produced different snapshots; docker's own output is
  non-deterministic. Investigate the tool/test; regeneration will not fix it.
- **`build_failed`** -- the Conda env or Singularity image failed to build before testing. Infra:
  build the env/image, then re-run (it blocks triage of that profile).
- **`no_ground_truth`** -- docker established no snapshot for the non-docker profiles to validate
  against (usually docker itself failed to produce one).
- **`syntax_error`** -- the Nextflow script failed to compile. Fix the `.nf`.
- Housekeeping statuses you may also see: `skipped`, `timeout` (exceeded the per-run timeout),
  `no_snapshot`, `n/a`.
- **`undeclared_outputs`** -- the tool produced files not declared in the module's `results`,
  `logs`, `versions`, or `nf_logs`. For each file help the user route it:
  - **`results`**: a real tool output users want (report, summary, data file)
  - **`logs`**: stdout/stderr from the tool
  - **`.outputs-ignore`**: staging artifact, intermediate, version-info side effect, or DB file
  `.outputs-ignore` lives at `modules/{name}/tests/.outputs-ignore` (one glob per line; `#`
  comments and blanks allowed; `staging/**` is ignored by default). NOTE: this check only runs
  when the tool **succeeds**, so it is masked on a profile that `tool_error`'d -- use
  `undeclared_outputs_union` to see the full set.
- **`tool_error`** -- the tool crashed at runtime. Read **`error_class`**:
  - `env_dependency` -- conda/singularity re-solved a too-new interpreter/dependency
    (e.g. py>=3.12 `pkg_resources`, numpy2 `newshape`, biopython `SeqFeature.strand`, R
    `readr`/`lifecycle` `deprecate_stop`). **Fix the env/recipe, NOT the test.**
  - `tool_crash` / `staging_bug` / `fs_permission` / `unknown` -- fix the module/upstream or
    the workspace; when `unknown`, read the `Command error:` block.
  `reason` carries the real tool error (from `Command error:`), not the downstream nf-test
  `NullPointerException`.

`generate` gates interpretation (shown in Run Parameters and the `# generate=<bool>` header
of `summary.tsv`):
- **`generate=true`**: the `.snap` was regenerated under docker first, so docker passing is the
  re-baseline; any drift shown is genuine (docker vs profile). `snapshot_stale` cannot occur.
- **`generate=false`**: docker **also** failing => `snapshot_stale` (run `--generate`). docker
  passing while a profile drifts => genuine content drift.

## Multi-profile layout & `summary.json`

Structured results live at `logs/run-tests/{ts}/summary.json` (plus `summary.tsv`, whose first
line is `# generate=<bool>`). Prefer `summary.json` for machine-readable detail; the CLI text is
the human summary. Key fields:

- `profiles[]`, `reference_profile` (`"docker"`).
- `results[].cells.{profile}`: `status`, `duration`, `reason`, `error_class` (tool_error only),
  `undeclared_outputs[]`.
- `results[].undeclared_outputs_union`: undeclared files unioned across profiles (unmasks
  profiles that `tool_error`'d).
- `results[].files[]`: per output file, the **cross-profile md5 matrix** -- `process`, `scope`
  (`sample`/`run`; subworkflow multi-record), `field`, `name`, `md5:{profile -> hash|null}`,
  `verdict`, `divergent_profiles[]`, plus:
  - `verdict`: `stable` (equal across all profiles that ran) | `divergent` | `indeterminate`
    (a profile didn't produce it) | `skip`.
  - `comparable`: `false` = intrinsically non-hashable (gz / normalized -> byte md5 is
    meaningless) => bucket **existence-only**; `verdict:"skip"`.
  - `incomplete[]`: profiles that produced no file (e.g. a `tool_error`'d conda) => `verdict`
    is `indeterminate`, NOT a false `stable`; re-check after fixing that profile.
  - `kind:"versions"` + `tool_key`: a `versions.yml` -> bucket to `contains('<tool_key>')`.
  This matrix is computed from actual runtime outputs, so it is populated **even on passing or
  stale cells** -- an always-on divergence diagnostic (also useful for `add-*` at creation time).
- `results[].suggested_edit` (module/subworkflow only): the exact test change implied by the
  verdicts -- `snapshot:[fields]` (stable), `existence:[fields]` (divergent content),
  `contains:[{field,value}]` (divergent versions). Subworkflow fields are **scope-prefixed**
  (`sample.`/`run.`, e.g. `sample.blast`, `run.versions`). Directly consumable and self-verifying
  (diff against the committed test). Workflow tier is intentionally `.nftignore`-only, so it has no
  `suggested_edit`; add the divergent globs to `workflows/**/tests/.nftignore` instead.

## Diagnostic files (per profile)

Layout: `logs/run-tests/{ts}/{tier}/{component}/{profile}/`:
- **`stdout.txt`** -- nf-test console, including the tool's own `Command error:` block.
  Read this for **`tool_error` root cause**.
- **`stderr.txt`** -- nf-test assertions, including the `Different Snapshot` per-file md5 diff.
  Read this for **drift / assertion detail**.
- **`outputs.txt`** -- `# Undeclared outputs:` list, or `# OK`.
- **`.nf-test/**`** -- preserved work tree (present for all cells, passing included), including
  `meta/output_0.json` (record field -> output file paths) and `meta/nextflow.log`.

Both `stdout.txt` and `stderr.txt` matter now, split by class (this replaces the old
"read stdout, not stderr" rule).

## Progressive Disclosure

Keep the initial summary compact and scannable. Do NOT open `stdout`/`stderr`/`nextflow.log`
during the initial summary -- the status matrix, `reason`, `error_class`, and `files[]` usually
suffice. When the user asks for deeper detail:

- **Specific component**: read its `summary.json` `results[]` entry first (cells, `reason`,
  `error_class`, `files[]`, `suggested_edit`). Then, if needed, open
  `{tier}/{component}/{profile}/stdout.txt` (tool_error) or that same dir's `stderr.txt` (drift diff).
- **Undeclared outputs**: use `undeclared_outputs_union` (or a cell's `undeclared_outputs[]`),
  then read the module's `main.nf` output block to advise `results` / `logs` / `.outputs-ignore`.
- **Tool / abort errors**: read the `Command error:` block in `{profile}/stdout.txt`; the full
  Nextflow log is at `{tier}/{component}/{profile}/.nf-test/tests/*/meta/nextflow.log`
  (focus on ERROR/WARN and the last ~50 lines).
- **Drift bucketing**: use `files[]` + `suggested_edit`; cross-check with the
  `Different Snapshot` block in `{profile}/stderr.txt`.

## Important Reminders

- CRITICAL: NEVER suggest `--update-snapshots` / snapshot regeneration for the drift statuses
  (`output_drift`/`version_drift`/`version+output_drift`) or env-drift `tool_error`s. Regen does NOT fix profile divergence --
  migrate the test (sccmec pattern) or fix the env. `--generate` is the fix **only** for
  `snapshot_stale`.
- `error_class: env_dependency` => fix the conda/singularity env or bioconda recipe, NOT the test.
- Read `{profile}/stdout.txt` for `tool_error` root cause and `{profile}/stderr.txt` for
  drift/assertion diffs -- both matter.
- `undeclared_outputs` can be masked on a `tool_error`'d profile -- always check
  `undeclared_outputs_union`.
- `files[]` verdicts: `comparable:false` (gz/normalized) -> existence-only; `verdict:indeterminate`
  + `incomplete:[...]` -> a profile didn't run, re-check after fixing (never a clean bill).
- Not all tiers/profiles appear in every run; a component with no Galaxy image has
  `galaxy:false` and no `singularity_galaxy` cell.
- `.nf-test/` work dirs are preserved per profile for **all** cells (passing included), so you can
  inspect any profile's `meta/output_0.json` or work tree -- not just failures.

## Updating Baselines

Baselines file: `conf/test-times.json`. Durations are **docker-based** (the CLI reports
"Docker duration").

To update baselines after a clean all-pass run, add `--update-baselines`:
```
bash .agents/skills/review-tests/scripts/run-bactopia-review-tests.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --silent --update-baselines
```
This writes actual runtimes from the current run into the baselines file and updates the
`_meta.updated` timestamp. Only entries for tested components are updated; other tiers
are left unchanged. After updating, re-run without `--update-baselines` to confirm anomalies
are resolved.

## Interpreting Timing Anomalies

Timing is measured against the **docker** profile.
- **generate=true vs generate=false**: a `generate=true` run executes tests twice (generate
  snapshots, then test against them). If baselines were recorded from a `generate=true` run but
  the current run uses `generate=false`, tests run at ~0.5x baseline -- expected, not suspicious.
- **Slow tests**: may reflect newly added test cases rather than regressions. Check recent
  commits to the component's test file before flagging.
- **Only flag anomalies as concerning** when the `generate` parameter matches between the
  baseline run and the current run.

## Self-Improvement

If you find yourself writing ad-hoc Python or bash to parse, explore, or extract data from the
CLI output or `summary.json`, that logic should be added to this skill or the underlying
`bactopia-review-tests` CLI instead. Update the skill so future sessions don't reinvent it.

## JSON Output

`logs/run-tests/{ts}/summary.json` is the primary structured source (schema above). The CLI can
also emit it with `--json` (add `--pretty` for readable output). See
`bactopia-review-tests --help` for details.
