---
name: review-tests
description: Review nf-test run results and present a diagnostic summary with grouped error analysis. Use when asked to review tests, check test results, show test failures, analyze test output, investigate why tests failed, see what's broken, or check test status. Accepts an optional timestamp argument to review a specific run.
---

# Review Tests

Analyze an nf-test run and present a diagnostic summary with grouped error analysis.

## Steps

### 1. Find the test run

If the user provided a timestamp argument (e.g., `/review-tests 20260324_081306`), look for `logs/{timestamp}/summary.json`. If it doesn't exist, list available runs from `logs/` and stop.

If no argument was provided, list directories in `logs/` at the repo root. Directory names are `YYYYMMDD_HHMMSS` timestamps -- take the lexicographically last one. If no directories exist, report "No test runs found in logs/" and stop.

### 2. Read summary.json

Read `logs/{timestamp}/summary.json`. Parse:
- `summary` object: status counts (passed, tool_error, non_reproducible, etc.)
- `results` array: per-component objects with component, tier, status, duration

### 3. Present the header

Display:
- Test run timestamp (reformat directory name to `YYYY-MM-DD HH:MM:SS UTC`)
- Total tests run (sum of all status counts)
- Tiers tested (check which tier subdirectories exist: modules, subworkflows)
- Duration: total (sum all durations), average per test, longest test (name + duration)

### 4. Status breakdown

Show overall counts as a compact table:

```
Status              Count
passed              70
tool_error          14
non_reproducible    1
```

### 5. Passed tests

Report only count and percentage: "70 passed (82.4% of 85 total)". Do NOT list individual passed components.

### 6. Investigate failures

For each non-passed status, read the `.stdout.txt` files from `logs/{timestamp}/{tier}/{component}.stdout.txt`. Process in priority order:

#### a. syntax_error

Compilation failures. Extract the error message from stdout. Show component name, tier, and the compilation error.

#### b. tool_error

Read each failed component's stdout file. Classify into error groups and present grouped:

**Group 1 -- Channel wiring errors**: Look for `java.lang.RuntimeException: Workflow has no output channels` or `ERROR ~ No such variable`. Extract the missing property and source file/line from `Check script '...' at line: N`. Present as: "Channel wiring errors (N): component1, component2 -- MissingPropertyException on ChannelOut"

**Group 2 -- Process execution failures**: Look for `ERROR ~ Error executing process > '...'`. Extract: process name, "Caused by:" line, exit status, and relevant lines from "Command error:" section. Present grouped by similar cause (e.g., download failures together, missing output files together).

**Group 3 -- Test assertion failures**: Look for `Assertion failed:` with `N of M assertions failed` but NO `ERROR ~` line. The workflow completed but test assertions failed. Present as: "Test assertion failures (N): component1, component2 -- workflow completed but output didn't match expectations"

**Group 4 -- Assertion with Nextflow error**: Both `Assertion failed:` AND `ERROR ~` present. Extract the ERROR line. Present grouped by error type.

#### c. non_reproducible

CRITICAL: NEVER suggest "rerun with --update-snapshots" -- that does NOT fix the root cause.

Report:
- Component name and tier
- The tool produced different output between two identical runs (generate mode ran the test twice)
- Suggest investigating WHY the output differs: non-deterministic tool behavior, timestamps in output, random seeds, floating-point ordering, thread-dependent output ordering
- Read the stdout file to see which assertions failed for additional context

#### d. snapshot_mismatch

Existing snapshot doesn't match current output. Could be expected if code was intentionally changed, or a regression. List components.

#### e. no_snapshot

Missing snapshot file. List components. Note these need a snapshot generated.

#### f. skipped

List components briefly.

## Error Extraction from stdout Files

The stdout files contain ANSI escape codes (`[31m`, `[0m`, etc.). Ignore these when parsing.

Key markers to search for:
- `ERROR ~ ` -- Primary Nextflow error message
- `Caused by:` -- One-line error reason
- `Assertion failed:` -- nf-test assertion failure count
- `Command exit status:` -- Process exit code
- `Command error:` -- Container stderr (useful for process failures)
- `Container:` -- Docker image used
- `Work dir:` -- Path for manual debugging
- `Check script '...' at line:` -- Source location of the error
- `Check '...nextflow.log'` -- Path to detailed Nextflow log

## Deeper Investigation

Do NOT read nextflow.log files during the initial summary. Only investigate when the user asks for deeper detail on a specific failure.

The path to the nextflow.log is in the stdout file after `Check '...' file for details`. These files only exist for failed tests (passed tests have `.nf-test/` cleaned up).

When reading nextflow.log (on request), focus on:
- Lines containing ERROR or WARN
- The last 50 lines (error/shutdown section)
- Stack traces indicating root cause

## Important Reminders

- Always read `.stdout.txt` files for diagnostics, NOT `.stderr.txt` -- stderr is misleading and lacks useful context
- The `logs/{timestamp}/` directory contains tier subdirectories based on what was tested -- not all tiers are present in every run
- Group errors by pattern for a scannable summary, not individual component dumps
- The `.nf-test/` work directories under component test dirs only exist for failed tests
- If the user asks about a specific component, offer to read its stdout file in full and check for nextflow.log
