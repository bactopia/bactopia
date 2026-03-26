---
name: review-tests
description: Review nf-test run results and present a diagnostic summary with grouped error analysis. Use when asked to review tests, check test results, show test failures, analyze test output, investigate why tests failed, see what's broken, or check test status. Accepts an optional timestamp argument to review a specific run.
---

# Review Tests

Run the review-tests CLI and interpret the JSON output for the user.

## Steps

1. Run `bactopia-review-tests` via the wrapper script:
   ```
   bash .claude/skills/review-tests/scripts/run-bactopia-review-tests.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --json --silent
   ```
   If the user provided a timestamp argument (e.g., `/review-tests 20260324_081306`),
   add `--run 20260324_081306`.

2. Parse the JSON output and present a clean, readable summary. Lead with:
   - Test run timestamp, total tests, tiers tested
   - Duration stats (total, average, longest test name + time)
   - Status breakdown as a compact table (counts per status)
   - Pass rate: "N passed (X% of M total)"

3. Present failure groups in priority order (skip groups with count 0):

   **undeclared_parameter**: Show the `parameters` summary (parameter name -> count).
   List all affected components grouped by parameter name for compactness.

   **missing_config**: Show each component with the `missing_path` from the message.

   **process_failure**: Show process name, `caused_by`, `exit_status`, and
   `command_error` excerpt. Group by similar cause if multiple.

   **null_container**: Note which container resolved to null.

   **abort_error**: Note nextflow.log path is available for deeper investigation.

   **assertion_failure**: Show `assertions_failed` / `assertions_total` per component.
   Note these likely need snapshot regeneration or investigation.

   **syntax_error**: Show component and compilation error message.

   **unclassified**: Show component and snippet for debugging.

4. If `timing_anomalies.baselines_available` is true:
   - Show slow tests (actual vs expected, ratio) if any
   - Show suspiciously fast tests if any
   - If none, note all tests within expected range

   If baselines are not available, mention `--update-baselines` can create them.

5. Summarize actionable items at the end:
   - Count of unique error patterns
   - Suggested next steps (e.g., "Declare `test_ont` parameter or remove from test params",
     "Add missing module.config files")

## Progressive Disclosure

The initial summary should be compact and scannable. When the user asks for deeper detail:

- **Specific component**: Read its stdout file at
  `logs/{timestamp}/{tier}/{component}.stdout.txt` using the Read tool
- **Abort errors**: Read the `nextflow_log` path from the JSON output
  (focus on ERROR/WARN lines and last 50 lines)
- **Assertion details**: Read the stdout file and look for specific assertion
  mismatch information

Do NOT read nextflow.log or stdout files during the initial summary.

## Important Reminders

- CRITICAL: NEVER suggest "rerun with --update-snapshots" for non_reproducible
  failures -- that does NOT fix the root cause
- Always read `.stdout.txt` files for diagnostics, NOT `.stderr.txt`
- The `logs/{timestamp}/` directory contains tier subdirectories based on what
  was tested -- not all tiers are present in every run
- The `.nf-test/` work directories under component test dirs only exist for
  failed tests
- If the user asks about a specific component, offer to read its stdout file
  in full and check for nextflow.log

## JSON Output Fields

- `timestamp` -- when the test run started (reformatted from directory name)
- `run_dir` -- path to the test run directory
- `total_tests` -- number of tests run
- `tiers_tested` -- list of tiers present (modules, subworkflows, workflows)
- `duration` -- `total_seconds`, `average_seconds`, `median_seconds`,
  `longest` (component/tier/duration), `shortest` (component/tier/duration)
- `status_counts` -- map of status to count (e.g., passed, tool_error)
- `passed` -- `count` and `percentage`
- `failure_groups` -- array of groups, each with:
  - `pattern` -- error classification key
  - `label` -- human-readable group name
  - `count` -- number of failures in this group
  - `detail` -- one-line explanation of this error type
  - `parameters` -- (undeclared_parameter only) map of param name to count
  - `components` -- array of affected components with pattern-specific fields:
    - Common: `component`, `tier`, `status`, `duration`, `pattern`, `message`
    - undeclared_parameter: `parameters` (list of param names)
    - missing_config: `missing_path`
    - process_failure: `process`, `caused_by`, `exit_status`, `command_error`
    - abort_error: `nextflow_log` (path for deeper investigation)
    - assertion_failure: `assertions_failed`, `assertions_total`
    - unclassified: `snippet` (first 500 chars)
- `timing_anomalies` -- `baselines_file`, `baselines_available`,
  `slow_tests` and `fast_tests` arrays with `component`, `tier`,
  `actual_seconds`, `expected_seconds`, `ratio`
