---
name: run-tests
description: Run Bactopia nf-tests via bactopia-test and produce a timestamped logs/ directory that /review-tests can interpret. Use when asked to run tests, execute tests, test a module, test a subworkflow, test a workflow, or validate a change. Accepts optional tier and component arguments (e.g., "run tests on snippy", "test abricate_run module", "test amrfinderplus subworkflow").
---

# Run Tests

Run the Bactopia nf-test suite through `bactopia-test` for a specific component
and present the live output to the user. This is the "before" half of the
`run-tests` / `review-tests` pair: this skill **runs** the tests and writes a
timestamped `logs/{timestamp}/` directory; `/review-tests` then **interprets**
that directory (grouping failures, reading stdout files, etc.). Keep the two
responsibilities clearly separated -- do not try to do `/review-tests`' job here.

## Steps

1. **Resolve `--tier` and `--include` from what the user said.** Use the
   "Resolving arguments" table below. If the user did not name a specific
   component, **stop and ask** which component they want to test -- do not run
   the full suite (see the guardrails block).

2. **Invoke the wrapper script** with the resolved arguments:

    ```
    bash .claude/skills/run-tests/scripts/run-bactopia-test.sh \
        --bactopia-path /home/rpetit3/repos/bactopia/bactopia \
        --test-data /home/rpetit3/repos/bactopia/bactopia-tests \
        --outdir /home/rpetit3/repos/bactopia/bactopia \
        --profile docker \
        --keep \
        [--tier TIER] \
        --include COMPONENT
    ```

    Omit `--tier` entirely when the user said "test snippy" (no tier hint) --
    the CLI's default `all` will then search every tier and the underscore-
    segment matcher (see "`--include` matching") will pick up every component
    whose name contains `snippy`.

3. **Present the CLI's live output directly.** `bactopia-test` produces a Rich
   table with per-component status, durations, and a final summary. Relay it
   without reformatting. Do not parse JSON, do not re-tabulate, do not read the
   stdout files the CLI writes.

4. **After the run finishes**, extract the run timestamp and hand off to
   `/review-tests`. See the "After the run" section.

## Resolving arguments (tier inference)

`bactopia-test`'s `--include` matcher is **exact-match OR underscore-segment
match** against a flattened component name (e.g.
`modules/snippy/run/tests/main.nf.test` → `snippy_run`). Patterns and wildcards
are **not** supported. That means `--include snippy` matches `snippy_run`,
`snippy_core`, etc.; `--include snippy_run` matches only the module.

Translate user phrasing into flags as follows:

| User says                                          | `--tier`       | `--include`    | Notes                                                                                                  |
| -------------------------------------------------- | -------------- | -------------- | ------------------------------------------------------------------------------------------------------ |
| "run all tests on snippy" / "test snippy"          | *(omit)*       | `snippy`       | Omitting `--tier` lets the CLI default (`all`) + segment match hit `snippy_run`, `snippy_core`, etc.   |
| "test abricate_run module" / "run tests on abricate_run" | `modules`      | `abricate_run` | `_run`/`_download`/`_predict` suffixes are a strong "this is a module" signal; pin `--tier modules`.   |
| "test amrfinderplus subworkflow"                   | `subworkflows` | `amrfinderplus`| The explicit word "subworkflow" pins the tier.                                                         |
| "test the snippy workflow"                         | `workflows`    | `snippy`       | The explicit word "workflow" pins the tier.                                                            |
| "run tests" / "run the test suite" / no component  | **REFUSE**     | —              | Stop and ask which component. The CLI default runs the full suite -- see guardrails.                  |
| "run all tests" (explicit)                         | **REFUSE, ASK**| —              | Ask the user which specific component; full-suite runs are slow and not what this skill is for.       |

If the user's request is ambiguous (e.g. "test snippy" could mean the module,
the subworkflow, or both), omitting `--tier` is the safe default -- the
segment matcher will catch all three and the user sees everything at once.

## Defaults the skill applies automatically

These flags are always added without asking the user:

| Flag              | Value                                        | Why                                                                           |
| ----------------- | -------------------------------------------- | ----------------------------------------------------------------------------- |
| `--bactopia-path` | `/home/rpetit3/repos/bactopia/bactopia`      | Canonical repo location on this machine.                                      |
| `--test-data`     | `/home/rpetit3/repos/bactopia/bactopia-tests`| Canonical test-data location; sets `BACTOPIA_TESTS`.                          |
| `--profile`       | `docker`                                     | Default execution profile. Docker is the baseline for reproducible tests.     |
| `--keep`          | *(always)*                                   | Preserves `.nf-test/` dirs and logs on pass; `/review-tests` needs them.      |
| `--outdir`        | `/home/rpetit3/repos/bactopia/bactopia`      | So `logs/{timestamp}/` lands at the repo root, where `/review-tests` reads.   |

## When to ask the user first (never auto-fill)

| Flag                 | Policy                                                                                                                                                                      |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--condadir`         | Only add if the user explicitly asked for `--profile conda`. **Ask the user** for the cache path before running; do not guess.                                             |
| `--singularity_cache`| Only add if the user explicitly asked for `--profile singularity`. **Ask the user** for the cache path before running; do not guess.                                       |
| `--generate`         | **Never** add unless the user's request contains one of: `--generate`, "generate mode", "regenerate snapshots", or "update snapshots". Destructively deletes snapshots.    |
| `--jobs N`           | Pass through if the user specified a number (e.g. "with 32 jobs"); otherwise omit and let the CLI default (`max(1, nproc // 4)`) apply.                                    |

## Important Reminders

These are the non-negotiable rules. Violating any of them can burn hours of
time or delete work the user cares about.

- **CRITICAL: never run with no component filter.** The CLI defaults to
  `--tier all` with no `--include`, which runs the entire test suite across
  96 modules + 87 subworkflows + 69 workflows. That is expensive and is not
  what this skill is for. If the user asks to "run tests" without naming a
  component, **stop and ask** which module / subworkflow / workflow they want.

- **NEVER pass `--generate`** unless the user explicitly asked for it with
  the literal flag name or the phrases "generate mode", "regenerate
  snapshots", or "update snapshots". `--generate` deletes existing snapshot
  files and re-runs each test twice. The deletion is not recoverable from
  within the skill.

- **ALWAYS pass `--keep`.** On pass, the CLI cleans up `.nf-test/` work
  directories by default; `--keep` preserves them. `/review-tests` may need
  to read into them to diagnose undeclared outputs or assertion mismatches.

- **ALWAYS pass `--outdir /home/rpetit3/repos/bactopia/bactopia`** so that
  `logs/{timestamp}/` is written at the bactopia repo root. `/review-tests`
  looks for logs relative to `--bactopia-path`; if `--outdir` is omitted the
  logs land in whatever directory the shell was invoked from and the
  downstream skill will not find them.

- **Do NOT pass `--json`.** The default Rich text output is human-friendly
  and already auto-logs a machine-readable `summary.json` into the logs
  directory. `--json` exists as a fallback for programmatic parsing -- it is
  not needed here.

- **Do NOT pass `--fail-fast`** unless the user explicitly asks. The default
  "run everything, report all failures at the end" behavior is what
  `/review-tests` expects to consume.

- **Do NOT interpret failures in detail here.** Do not read
  `.stdout.txt`, `.stderr.txt`, `.outputs.txt`, or `nextflow.log` files. Do
  not group failures by type. Do not recommend fixes. That is deliberately
  reserved for `/review-tests` so the two skills stay loosely coupled and
  each has a single clear job.

- **Profile switching requires extra info.** If the user asks for
  `--profile singularity` or `--profile conda`, switch profiles **and** ask
  them for the corresponding `--singularity_cache` or `--condadir` path
  before invoking the wrapper. Do not guess.

## After the run

When `bactopia-test` finishes, do these four things -- nothing more:

1. **Report overall pass/fail counts** from the CLI's final summary table
   (it prints a row-count breakdown by status).

2. **Extract the run timestamp.** The CLI prints the path to the logs
   directory, which ends in a `YYYYMMDD_HHMMSS` directory (e.g.
   `logs/20260410_143022/`). Pull that timestamp out and show it to the user.

3. **Point the user at `/review-tests`** with an exact next step:

    > Run `/review-tests {timestamp}` for a diagnostic grouping of any
    > failures and per-component details.

4. **Stop.** Do not open log files. Do not enumerate failed components
   beyond the CLI's own summary table. Do not re-run. The job of this skill
   is to produce the logs directory and hand off.

## Notes

### `--include` matching (full semantics)

`bactopia-test` builds a flattened component name by taking the test file
path, stripping the `tests/main.nf.test` suffix and the tier prefix, and
joining the remaining segments with `_`. Examples:

- `modules/snippy/run/tests/main.nf.test` → `snippy_run`
- `modules/bactopia/gather/tests/main.nf.test` → `bactopia_gather`
- `subworkflows/snippy/core/tests/main.nf.test` → `snippy_core`
- `subworkflows/bactopia-tools/amrfinderplus/tests/main.nf.test` → `amrfinderplus`
  (the `bactopia-tools/` prefix is stripped)
- `workflows/teton/tests/main.nf.test` → `teton`

`--include X` matches a component if:

- `X` equals the full flattened name exactly, OR
- `X` appears as a segment when the flattened name is split on `_`.

Concrete consequences:

- `--include snippy` matches `snippy_run`, `snippy_core`, and `snippy`
  (any tier) but **not** `snippyassembly` -- no underscore boundary.
- `--include run` would match every `*_run` component, which is almost
  never what the user wants. Prefer the more specific form.
- Wildcards, globs, and regex are not supported. Multiple values are
  comma-separated (`--include snippy,bakta`).

### CLI flag reference (complete list)

Taken from [bactopia-py/bactopia/cli/testing.py](/home/rpetit3/repos/bactopia/bactopia-py/bactopia/cli/testing.py).
Defaults in parentheses.

**Paths (required)**
- `--bactopia-path` — bactopia repo root (required)
- `--test-data` — test-data directory, sets `BACTOPIA_TESTS` (required unless `--cleanup`)

**Selection**
- `--tier` — `modules` / `subworkflows` / `workflows` / `all` (default: `all`)
- `--include` — comma-separated component names (default: none → all)
- `--exclude` — comma-separated component names (default: none)

**Execution**
- `--profile` — `docker` / `singularity` / `conda` (default: `docker`)
- `--condadir` — conda cache dir (default: `${BACTOPIA_CACHEDIR}/conda`)
- `--singularity_cache` — singularity cache dir (default: `${BACTOPIA_CACHEDIR}/singularity`)
- `--generate` — delete snapshots and run twice (default: off)
- `--jobs` — parallel workers (default: `max(1, nproc // 4)`)
- `--fail-fast` — stop on first failure (default: off)
- `--timeout` — per-test timeout in **minutes**, 0 to disable (default: 90)

**Cleanup (operates instead of running tests)**
- `--cleanup` — find and remove all `.nf-test/` temp files, then exit
- `--dry-run` — with `--cleanup`, list what would be removed

**Output**
- `--outdir` — directory to write `logs/` into (default: `.`)
- `--keep` — preserve `.nf-test/` dirs and logs on pass (default: off)
- `--json` — emit results as JSON (default: off; used only as fallback)

**Logging**
- `--verbose` — DEBUG logging
- `--silent` — ERROR logging only
- `--version` / `--help`

### Output layout written by the CLI

```
{outdir}/logs/{YYYYMMDD_HHMMSS}/
├── summary.json                               # machine-readable rollup
├── summary.tsv                                # same data in TSV
├── modules/
│   ├── {component}.stdout.txt                 # captured stdout from nf-test
│   ├── {component}.stderr.txt                 # captured stderr
│   └── {component}.outputs.txt                # undeclared-outputs report
├── subworkflows/
│   └── ... (same structure)
└── workflows/
    └── ... (same structure)
```

Only the tiers that were tested have subdirectories in a given run.
`/review-tests` reads these files directly -- do not pre-load them here.

### Wrapper script discovery order

The wrapper at `.claude/skills/run-tests/scripts/run-bactopia-test.sh`
locates `bactopia-test` by checking, in order:

1. `bactopia-test` on `PATH` (respects an already-activated env)
2. A conda env whose path ends in `/bactopia-dev` (exact) -- **preferred**,
   built from `environment.yml` at the bactopia repo root
3. A conda env whose path ends in `/bactopia-py` (exact) -- stable release
4. A conda env whose path contains `/bactopia-` (fuzzy, first match) --
   catches `bactopia-py-dev`, `bactopia-dev-v4`, and any other
   `bactopia-*` env
5. Fails with an error directing the user to install `bactopia-py` or
   activate an env containing it.

This ordering matters: `bactopia-dev` is the canonical development env
(nextflow, nf-test, bactopia-py, and nf-bactopia build tools all in one
place) and will usually have the freshest `bactopia-test` version. The
older `bactopia-py-dev` env -- which contained only the Python CLI -- is
still supported as a fallback via the fuzzy step but should not be
preferred over `bactopia-dev`.

If a conda env is selected, the wrapper uses `conda run --prefix <path>`
to invoke the CLI without requiring the env to be activated. All arguments
are forwarded through `"$@"`.

### Sibling skills

- `/review-tests` — the "after" half. Reads `logs/{timestamp}/`, groups
  failures by type, reads stdout files on request, and suggests next steps.
  Always point the user here after a run completes.
- `/project-status` — component counts and coverage. Unrelated to the test
  loop but uses the same wrapper-script pattern.
- `/update-module` — bumps tool versions. Unrelated, but is the reference
  for the "ask the user before mutating" pattern borrowed here for
  `--generate` and profile switches.
