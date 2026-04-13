---
name: review-docs
description: Review staleness of reference docs under .claude/docs/ using bactopia-docs --validate. Detects deprecated patterns (residue from past migrations like flattenPaths, the 4-channel emission framing, meta:Map) and ground-truth violations (stale module/subworkflow/workflow counts, wrong Nextflow version, references to nonexistent bactopia-* commands or lint rule IDs, broken markdown link targets). Use this skill whenever the user asks to review docs, check doc staleness, audit reference docs, find outdated documentation, verify doc claims, check if docs are current, or scan .claude/docs for drift after a migration.
---

# Review Docs

Run `bactopia-docs --validate` and present the staleness report, then triage findings on request. Cross-cutting (per-repo) doc audit, parallel to `/review-citations` for citations and `/review-groovydoc` for component-level GroovyDoc.

## What this covers

Two check families that don't fit the per-component lint rule model:

- **Deprecated patterns (D0xx)** — regex matches against the `data/docs-patterns.yml` registry. Each entry flags a phrase retired by a past migration (e.g. `flattenPaths`, the `4-channel` subworkflow emission framing, `Tuple<Map, ...>` pre-Record typing, `meta: Map`, `EMPTY_*` placeholders, `data/catalog.json` path, `bactopiatool_init` name). New migrations append entries to the YAML; the rule list grows over time.
- **Ground-truth assertions (D1xx)** — claims about repo state that are derivable from the live tree and SHOULD match what the docs say:
    - **D101/D102/D103** — module/subworkflow/workflow counts (`find <tier> -name main.nf | wc -l`)
    - **D104** — Nextflow version (`nextflowVersion` in `nextflow.config`); informational mentions like `26.04+` or `until Nextflow X` are skipped.
    - **D105** — `` `bactopia-*` `` references inside backticks must resolve to a `[tool.poetry.scripts]` entry in `bactopia-py/pyproject.toml`. Bare prose mentions (`bactopia-tools`, `bactopia-py`) are ignored.
    - **D106** — `M0xx`/`S0xx`/`W0xx`/`MC0xx`/`JS0xx`/`FMT0xx` lint rule IDs must resolve to a `rid = "..."` assignment in `bactopia-py/bactopia/lint/rules/`.
    - **D108** — markdown link targets `[text](path)` must resolve to a real file. URLs and anchor-only links are skipped.

Component-level checks live elsewhere:
- Module/subworkflow GroovyDoc → `/review-groovydoc` (M0xx/S0xx rules, including `@citation` keys via M035/S019)
- Citation integrity (`citations.yml` orphans, workflow `@citation` references) → `/review-citations`

## Steps

1. Run `bactopia-docs --validate` via the wrapper, asking for JSON so it's easy to parse:

    ```
    bash .claude/skills/review-docs/scripts/run-bactopia-docs.sh \
        --bactopia-path /home/rpetit3/repos/bactopia/bactopia \
        --validate --json --silent
    ```

    Exit code is `0` when all checks pass, `1` when any FAIL is found. Both cases print JSON on stdout.

2. Parse the JSON. Shape:

    ```json
    {
      "bactopia_path": "/home/rpetit3/repos/bactopia/bactopia",
      "docs_path": ".claude/docs",
      "patterns_file": "data/docs-patterns.yml",
      "ground_truth": {
        "counts": {"modules": 97, "subworkflows": 88, "workflows": 70},
        "nextflow_version": "25.04.6",
        "cli_commands_total": 28,
        "lint_rule_ids_total": 104,
        "bactopia_py_resolved": "/home/rpetit3/repos/bactopia/bactopia-py"
      },
      "files_scanned": ["standards/01-style-guide.md", "..."],
      "deprecated_patterns": [
        {
          "rule_id": "D001", "severity": "FAIL",
          "file": "reference/03-glossary.md", "line": 40,
          "match": "**flattenPaths**",
          "pattern": "flattenPaths",
          "hint": "Function removed from nf-bactopia plugin; remove or rephrase."
        }
      ],
      "ground_truth_violations": [
        {
          "rule_id": "D101", "severity": "FAIL",
          "file": "project/01-repository-structure.md", "line": 42,
          "match": "**Count**: 96 modules",
          "claim": "96 modules", "actual": "97 modules"
        }
      ],
      "summary": {
        "files_scanned": 16,
        "deprecated_pattern_hits": 11,
        "ground_truth_violations": 5,
        "fail": 16, "warn": 0,
        "patterns_loaded": 8
      }
    }
    ```

    Key fields:
    - `deprecated_patterns[]` — D0xx hits, one per (file, line, pattern) match.
    - `ground_truth_violations[]` — D1xx hits; carries either `claim`/`actual` (counts/version) or `reference`/`hint` (CLI/rule/path).
    - `ground_truth.bactopia_py_resolved` — `null` if the sibling repo wasn't found; D105/D106 are skipped silently in that case.

3. Present a summary grouped by check family:

    - **Clean repo** (`fail == 0`): report `"All <files_scanned> docs clean — <patterns_loaded> deprecated patterns, <cli_commands_total> CLI commands, <lint_rule_ids_total> lint rule IDs checked."`. Stop there.
    - **Issues found**: show sections in this order:
        1. **Deprecated patterns (D0xx)** — group by `rule_id` so the user sees migration residue together. For each rule, list `file:line — match` with the hint shown once at the top of the group.

            ```
            D001 (flattenPaths): function removed from nf-bactopia plugin
              • reference/03-glossary.md:40   - **flattenPaths**
              • reference/03-glossary.md:155  - **flattenPaths**: Convert Set<Path> to Path
              • reference/04-plugin-functions.md:114  ### flattenPaths (deprecated)
              ...
            ```

        2. **Ground-truth violations (D1xx)** — group by `rule_id`. For count/version checks, show `claim → actual`; for reference checks, show the bad reference + suggested fix.

            ```
            D101 (module count): claim 96 → actual 97
              • project/01-repository-structure.md:42

            D104 (Nextflow version): claim 26.01.0 → actual 25.04.6
              • project/04-testing-framework.md:4
              • project/04-testing-framework.md:316
            ```

        Lead with whichever family has more items. Keep it scannable: if any rule has more than ~10 hits, show the first 5 plus `(... N more)` and offer to print the rest on request.

4. If the user asks to fix issues, walk through them one rule at a time. Different rules need different treatments — don't batch them.

### Fixing deprecated patterns (D0xx)

Each pattern represents a migration that already landed; the doc is just out of sync. Three fix patterns:

- **Replace with current term**: most common case. `flattenPaths` → describe the current direct-emit pattern. `4-channel` → `2 channels (sample_outputs + run_outputs)`. `Tuple<Map, Path>` → `Channel<Record>`. Read the surrounding paragraph before editing — sometimes the whole sentence needs rewriting, not just a token swap. Look at [reference/01-examples.md](.claude/docs/reference/01-examples.md) for the modern equivalents.
- **Delete entirely**: glossary entries for terms that no longer exist (e.g. the `flattenPaths` definition, the `Tuple<Map, ...>` type entries) should be removed, not rephrased. Take the whole bullet/section.
- **Mark as historical with inline ignore**: very rarely, a doc legitimately needs to mention a deprecated term — e.g. a "what changed in v4" note. Suppress the rule on that line with an HTML comment:

    ```markdown
    Pre-v4 subworkflows used the 4-channel emission pattern. <!-- bactopia-docs: ignore D002 -->
    ```

    Multiple rules: `<!-- bactopia-docs: ignore D001, D002 -->`. Use this sparingly — every suppression is technical debt that will outlive the reason it was added.

### Fixing ground-truth violations (D1xx)

- **D101/D102/D103 (counts)**: simple integer swap. The CLI's `actual` value is authoritative — re-run after editing to confirm.
- **D104 (Nextflow version)**: replace the claimed version with the actual `nextflowVersion` from `nextflow.config`. If the doc was making a forward-looking claim ("targeting Nextflow X"), reword to remove the version or add `+`/`x` so D104 treats it as informational.
- **D105 (CLI references)**: either fix a typo (e.g. `bactopia-statuz` → `bactopia-status`) or remove the reference if the command genuinely doesn't exist. Don't invent commands.
- **D106 (lint rule IDs)**: same — fix typo or remove. The valid set lives in [bactopia-py/bactopia/lint/rules/](../../../bactopia-py/bactopia/lint/rules/).
- **D108 (broken links)**: either fix the path (compare against the actual repo layout) or remove the link if the target was deleted. URLs aren't checked, so external links never trip this.

After edits, re-run step 1. A clean re-run proves the fix landed correctly.

### Adding a new deprecated pattern

When a migration retires a term/pattern that should never appear again in docs, append it to [data/docs-patterns.yml](data/docs-patterns.yml):

```yaml
- id: D009                         # next available D0xx
  pattern: oldTermName             # regex by default
  literal: true                    # set if pattern should be matched as plain text
  severity: FAIL                   # PASS | WARN | FAIL
  hint: "What to do instead."
  rationale: |
    Free text describing the migration that retired this term.
```

Always confirm with the user before adding a pattern — false positives at the registry level multiply across every doc forever.

## Important constraints

- **Confirm before editing docs.** The CLI reports drift objectively, but choosing between "rewrite the section", "delete the entry", and "mark historical with inline ignore" is a judgment call. Pause and ask.
- **Don't fabricate facts.** If a count/version/reference is broken and the user doesn't know the right replacement, leave the doc alone and surface the gap. The CLI's `actual` field is reliable for counts/version; for CLI/rule references, the canonical source is the bactopia-py codebase.
- **D108 path resolution is heuristic.** It checks markdown link targets relative to the doc's directory, then falls back to the repo root. Symlinks and case-insensitive filesystems can produce edge-case false positives — always inspect the link before deleting it.
- **Suppression is a last resort.** `<!-- bactopia-docs: ignore Dxxx -->` is technical debt; prefer rewriting to remove the deprecated term entirely. The rare legitimate use is a deliberate historical reference (changelog, "what changed" notes).
- **D105/D106 skip silently when bactopia-py isn't found.** If `ground_truth.bactopia_py_resolved` is `null`, those checks didn't run. Pass `--bactopia-py-path` explicitly if the sibling repo is in a non-default location.
- The wrapper script auto-discovers `bactopia-docs` (checks PATH, then conda envs). No need to activate an env first.
- The CLI's `--bactopia-path` must point at the repo root so the validator can find `.claude/docs/`, `data/docs-patterns.yml`, `nextflow.config`, and the tier directories.

## Quick reference

### Inline suppression syntax

```markdown
Some line that legitimately mentions flattenPaths. <!-- bactopia-docs: ignore D001 -->
Multiple rules on one line. <!-- bactopia-docs: ignore D001, D002 -->
```

The HTML comment must be on the same line as the match. Suppression is per-rule, per-line — applies only to that line.

### `data/docs-patterns.yml` entry shape

```yaml
patterns:
  - id: D001                       # D001+, sequential
    pattern: flattenPaths          # regex by default; set literal:true for plain text
    literal: true                  # optional
    severity: FAIL                 # PASS | WARN | FAIL (default FAIL)
    hint: "Short remediation."     # shown in CLI output
    rationale: |                   # free-text history
      The migration / decision that retired this term.
```

### CLI flags

- `--bactopia-path PATH` — required, points at the bactopia repo root
- `--docs-path PATH` — relative to bactopia-path (default: `.claude/docs`)
- `--patterns-file PATH` — relative to bactopia-path (default: `data/docs-patterns.yml`)
- `--bactopia-py-path PATH` — override the sibling-repo discovery for D105/D106
- `--skip-path-check` — skip D108 (faster runs; useful when you know link health is fine)
- `--json` — emit structured output instead of rich tables
- `--silent` — suppress non-error output when the run is clean
- `--plain-text` / `-p` — disable rich formatting for piping

### When to redirect to other skills

- Module / subworkflow GroovyDoc accuracy → `/review-groovydoc`
- Citation integrity (`citations.yml` orphans, workflow `@citation` refs) → `/review-citations`
- Component coverage / project state → `/project-status`
- Catalog regeneration → `/update-catalog`
