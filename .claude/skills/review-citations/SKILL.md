---
name: review-citations
description: Review citation integrity across data/citations.yml and @citation tags using bactopia-citations --validate. Detects orphan citation keys (defined in the yml but never referenced) and workflow @citation keys that don't resolve to a yml entry. Use this skill whenever the user asks to review citations, check citation integrity, audit citations.yml, find orphan citations, clean up unused citations, validate workflow @citation tags, or verify that every tool cited in a workflow has a matching entry in the citations file.
---

# Review Citations

Run `bactopia-citations --validate` and present the integrity report, then apply fixes on request.

## What this covers

Three gaps that the per-component lint rules don't catch:

- **Orphan keys**: entries defined in `data/citations.yml` that are never referenced by any `@citation` tag. Accumulate as tools are replaced or renamed. The CLI suggests candidate "homes" for each orphan based on repo state — use those to decide between wire-up vs. drop.
- **Expected orphans**: yml entries marked `provenance_only: true` — foundational tools acknowledged for provenance only (e.g. `nextflow`, `nf-core`, `nf-test`). Surface them as informational, never flag as an issue.
- **Workflow missing keys**: `@citation` keys in workflow `main.nf` files that don't resolve to an entry in `citations.yml`. The W-series lint rules have no M035/S019 equivalent, so these only surface here.

Module and subworkflow `@citation` validation is handled by `bactopia-lint` rules M035 (modules) and S019 (subworkflows) — when the user wants to check those, point them at `/review-groovydoc`.

## Steps

1. Run `bactopia-citations --validate` via the wrapper, asking for JSON so it's easy to parse:

    ```
    bash .claude/skills/review-citations/scripts/run-bactopia-citations.sh \
        --bactopia-path /home/rpetit3/repos/bactopia/bactopia \
        --validate --json --silent
    ```

    Exit code is `0` when there are no real orphans and no missing workflow keys (expected orphans alone do not trigger a non-zero exit). Exit `1` when any real issue is found. Both cases print JSON on stdout.

2. Parse the JSON. Shape:

    ```json
    {
      "orphans": {
        "workflows": [],
        "datasets_ariba": ["srst2"],
        "tools": ["aragorn", "barrnap"]
      },
      "expected_orphans": {
        "influences": ["nfcore"],
        "tools": ["nextflow", "nftest"]
      },
      "potential_homes": {
        "aragorn": [
          {"type": "toolName", "path": "modules/prokka/module.config:30", "hint": "ext.toolName = \"bioconda::prokka=1.15.6\""},
          {"type": "script_token", "path": "modules/prokka/main.nf:55", "hint": "prokka --outdir ..."}
        ],
        "srst2": [
          {"type": "sibling_key", "path": "data/citations.yml (arg_annot)", "hint": "possible duplicate/variant of referenced key 'arg_annot'"}
        ]
      },
      "missing_workflow_keys": [
        {"component": "workflows/foo", "file": "workflows/foo/main.nf", "line": 42, "key": "typoed_key"}
      ],
      "summary": {
        "orphans_total": 3,
        "expected_orphans_total": 3,
        "missing_total": 1,
        "yml_total": 157,
        "referenced_total": 154
      }
    }
    ```

    Key fields:

    - `orphans` — real issues that need user judgment.
    - `expected_orphans` — `provenance_only: true` entries; informational only, never a failure.
    - `potential_homes[orphan]` — heuristic candidates for where each real orphan could be wired. Candidate types: `directory`, `toolName`, `config_param`, `script_token`, `sibling_key` (ranked roughly by specificity in that order).
    - `missing_workflow_keys` — workflow-tier `@citation` references that don't resolve.

3. Present a summary grouped by issue type:

    - **Clean repo** (orphans + missing both zero): report `"All <yml_total> citations are referenced, <expected_orphans_total> provenance-only acknowledged"` and list the expected-orphan keys inline as a dim note. Stop there.
    - **Issues found**: show sections in this order:
        1. **Orphans** — per orphan, render the key with its top 2-3 `potential_homes` candidates inline. Lead with the most specific candidate (directory > toolName > config_param > script_token > sibling_key). Example:

            ```
            aragorn (tools) — candidates:
              • toolName: modules/prokka/module.config:30  ext.toolName = "bioconda::prokka=1.15.6"
              • script_token: modules/prokka/main.nf:55  prokka --outdir ...
            ```

        2. **Missing workflow keys** — component, file:line, bad key.
        3. **Expected orphans** — single dim line listing the keys so the user knows they exist but aren't flagged.

        Lead with whichever of (1) or (2) has more items. Keep the output scannable; if there are many orphans, show counts + top 5-10 with homes, then note the rest by count.

4. If the user asks to fix issues, walk through them one category at a time. Orphans and missing keys need different treatment, so don't batch them together.

### Fixing orphans (three sub-cases)

Orphans need per-key judgment — don't sweep them. For each real orphan, decide between three actions:

- **Wire up**: The orphan has a plausible home in `potential_homes`. Edit the target `main.nf` to append the orphan key to its `@citation` line (keep existing cites first, append new keys). Verify the result with the Edit tool — if the current `@citation` line is `* @citation prokka`, the updated line becomes `* @citation prokka, aragorn`.
- **Drop**: Truly retired — no home, no sibling, no future plans. Use the Edit tool to delete the key line plus its entire indented block in `data/citations.yml` (typically 4-6 lines; `cite: |` multi-line bodies run longer).
- **Keep with provenance marker**: The orphan is deliberately unreferenced (foundational tool, planned but not yet wired, retained for historical record). Add `provenance_only: true` to its yml block — this moves it into `expected_orphans` on the next run so it stops tripping the validator. Place the field right below `name:`. Example:

    ```yaml
    nextflow:
      name: Nextflow
      provenance_only: true
      link: ...
      cite: |
        ...
    ```

Every `citations.yml` entry shape:

```yaml
aragorn:
  name: Aragorn
  link: ...
  description: ...
  cite: |
    Citation text across one or more lines.
```

Keys live under six top-level sections: `workflows`, `datasets_ariba`, `datasets_generic`, `datasets_minmer`, `influences`, `tools`. The `orphans` field in the JSON is already grouped by section — use that grouping so you're editing under the right parent key.

**Always confirm before editing.** Heuristics produce candidates, not decisions. Present the top suggestion and ask the user to pick wire-up / drop / keep-with-marker.

### Fixing missing workflow keys

Two sub-cases, determined by inspecting the bad key:

- **Typo** — the key looks like a variant of an existing key (e.g. `clonalframeml` vs `clonal_frame_ml`). Suggest the closest valid key from `citations.yml`. Edit the workflow `main.nf` to correct the `@citation` line. The file and line number are in the JSON.
- **Genuinely new tool** — the key is a real tool that just hasn't been added to the yml yet. Ask the user for the entry details (name, cite, optionally link/description) and add a new block in the appropriate section of `data/citations.yml`. Don't invent citations — if the user doesn't have the details at hand, leave the `@citation` as-is and flag it for follow-up.

After edits, re-run step 1 to confirm the issue is resolved. A clean re-run proves the fix landed correctly.

## Important constraints

- Always confirm before editing `citations.yml`. The CLI reports orphans objectively, but deciding between wire-up / drop / provenance-marker is a judgment call that needs the user's input.
- Don't fabricate citations. If a workflow references a tool that genuinely isn't in the yml and the user doesn't have the citation text, the right answer is to pause and flag it, not to paper over the gap.
- `potential_homes` is heuristic. A top candidate is a strong hint but not proof — always check the target `main.nf` before writing to it. Occasionally the "script_token" heuristic flags incidental string matches rather than real usage.
- Expected orphans are informational only. Don't pull them into the orphan triage loop unless the user explicitly asks about them (they already made the provenance decision when the marker was added).
- The wrapper script auto-discovers `bactopia-citations` (checks PATH, then conda envs). No need to activate an env first.
- The CLI's `--bactopia-path` must point at the repo root (e.g. `/home/rpetit3/repos/bactopia/bactopia`) so the walker can find `modules/`, `subworkflows/`, `workflows/`, and `data/citations.yml`.

## Quick reference

### `@citation` tag format

```
 * @citation key1, key2, key3
```

- Comma-separated, case-insensitive keys (stored lowercase in the yml)
- Must match a top-level key under one of the six sections in `data/citations.yml`
- Multi-line `@citation` continuations are supported by the GroovyDoc parser; validation treats them the same as a single line

### `citations.yml` sections

- `workflows` — top-level pipelines (bactopia, staphopia, ...)
- `datasets_ariba`, `datasets_generic`, `datasets_minmer` — reference datasets
- `influences` — foundational tools cited for inspiration (currently just nf-core)
- `tools` — most individual bioinformatics tools (150+ entries)

### `provenance_only: true` marker

Entries with this flag are intentionally unreferenced and surface in `expected_orphans` instead of `orphans`. Use it for foundational infrastructure tools that belong on the docs/citations page but don't have a natural home on any individual pipeline component. Current entries: `nfcore`, `nextflow`, `nftest`.

### When to redirect to other skills

- Module or subworkflow `@citation` key typos → `/review-groovydoc` (rules M035, S019)
- Workflow GroovyDoc structure issues unrelated to citations → no skill yet; manual review against `.claude/docs/standards/06-workflow-documentation.md`
