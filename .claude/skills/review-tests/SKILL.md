---
name: review-tests
description: Review nf-test run results and present a diagnostic summary with grouped error analysis. Use when asked to review tests, check test results, show test failures, analyze test output, investigate why tests failed, see what's broken, or check test status. Accepts an optional timestamp argument to review a specific run.
---

# Review Tests

Run the review-tests CLI and present the results to the user.

## Steps

1. Run `bactopia-review-tests` via the wrapper script using the **default text output**
   (do NOT use `--json`):
   ```
   bash .claude/skills/review-tests/scripts/run-bactopia-review-tests.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --silent
   ```
   If the user provided a timestamp argument (e.g., `/review-tests 20260324_081306`),
   add `--run 20260324_081306`.

2. Present the text output directly to the user. The CLI already produces a clean,
   well-formatted summary with tables. Do NOT parse JSON or write extra code to
   reformat -- just relay the output with your interpretation.

3. Add interpretation and context after showing the output:
   - For **assertion_failure** results, check the run parameters shown in the output:
     - If `generate` was **true**: snapshots were regenerated and the test was run a
       second time against them. These are **real failures** -- the workflow output
       does not match its own freshly-generated snapshot, meaning the output is
       non-deterministic or the test assertions are wrong. Flag these as needing
       investigation, NOT snapshot regeneration.
     - If `generate` was **false** (or not shown): snapshots may be stale.
       Note these likely need snapshot regeneration or investigation.
   - For **suspiciously fast tests**: note these likely exited early without running.
   - Summarize actionable items and suggested next steps.

4. If the text output is too large for a single response, summarize the key sections
   (overview, status breakdown, failures) and note that timing details are available
   on request. Use `--json` only as a fallback if the text output cannot be displayed.

## Progressive Disclosure

The initial summary should be compact and scannable. When the user asks for deeper detail:

- **Specific component**: Read its stdout file at
  `logs/{timestamp}/{tier}/{component}.stdout.txt` using the Read tool
- **Abort errors**: Read the nextflow.log for the component
  (focus on ERROR/WARN lines and last 50 lines).
  To find the log path, re-run with `--json` and check the `nextflow_log` field,
  or look in `logs/{timestamp}/{tier}/{component}.stdout.txt` for the path.
- **Assertion details**: Read the stdout file and look for specific assertion
  mismatch information

Do NOT read nextflow.log or stdout files during the initial summary.

## Important Reminders

- CRITICAL: NEVER suggest "rerun with --update-snapshots" for non_reproducible
  failures -- that does NOT fix the root cause
- When `params.generate` is true, NEVER suggest snapshot regeneration for
  assertion failures -- snapshots were already regenerated during this run.
  These represent non-deterministic output or incorrect test assertions.
- Always read `.stdout.txt` files for diagnostics, NOT `.stderr.txt`
- The `logs/{timestamp}/` directory contains tier subdirectories based on what
  was tested -- not all tiers are present in every run
- The `.nf-test/` work directories under component test dirs only exist for
  failed tests
- If the user asks about a specific component, offer to read its stdout file
  in full and check for nextflow.log

## JSON Output (Fallback)

The `--json` flag is available as a fallback for programmatic access or when the
text output is too large. Use it with `--pretty` for readable JSON. See
`bactopia-review-tests --help` for details on JSON fields.
