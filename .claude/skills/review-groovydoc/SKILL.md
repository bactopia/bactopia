---
name: review-groovydoc
description: Review GroovyDoc accuracy across modules and subworkflows using bactopia-lint. Checks @output/@input field matching, @modules/@subworkflows lists, citation keys, tag ordering, and formatting. Use when asked to review GroovyDoc, check documentation accuracy, validate module/subworkflow docs, or audit GroovyDoc.
---

# Review GroovyDoc

Run bactopia-lint focused on GroovyDoc accuracy rules and present the results.

## Steps

1. Run `bactopia-lint` via the wrapper script. By default, lint both modules and subworkflows:
   ```
   bash .claude/skills/review-groovydoc/scripts/run-bactopia-lint.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --quiet --json --silent
   ```
   - For modules only: add `--no-subworkflows --no-workflows`
   - For subworkflows only: add `--no-modules --no-workflows`
   - For a specific module: add `--module <name>`

2. Parse the JSON output. Filter results to **GroovyDoc-related rules only**:

   **Module rules (M-series):**
   - M006: GroovyDoc block present
   - M007: Required GroovyDoc tags
   - M008: @status value valid
   - M009: @tags sub-keys present
   - M010: Features list formatting
   - M014: Links use HTTPS
   - M031: @output fields match actual record() output
   - M032: @input record fields match actual input
   - M033: Optionality markers (?) match between GroovyDoc and code
   - M034: Standard fields not described in @output
   - M035: @citation keys exist in citations.yml
   - M036: GroovyDoc tag ordering
   - M037: Blank lines between GroovyDoc sections
   - M038: No `*/` inside GroovyDoc block

   **Subworkflow rules (S-series):**
   - S003: GroovyDoc with required tags
   - S004: Features comma-separated without spaces
   - S006: Links use HTTPS
   - S017: @modules match actual module includes
   - S018: @subworkflows match actual subworkflow includes
   - S019: @citation keys exist in citations.yml
   - S020: @tags complexity value is valid
   - S021: @tags input-type value is valid
   - S022: @tags output-type value is valid
   - S023: @tags features values are valid
   - S024: GroovyDoc tag ordering
   - S026: All emit channels have a corresponding @output tag
   - S027: @output field descriptions must not exist for channel.empty() emits

   **Entry workflows (W-series):**
   - No dedicated GroovyDoc-accuracy W-rules exist. Workflow GroovyDoc review is a manual pass — check the workflow main.nf against `.claude/docs/standards/06-workflow-documentation.md`, verify `@subworkflows` directory keys match actual `include { ... }` imports (e.g., `utils_bactopia-tools`, `bactopia_qc`), validate `@citation` keys against `data/citations.yml`, confirm `@input` names match declared `params { }`, and ensure each `@section` has ≥1 `@publish`. The one-off drift-check approach used on 2026-04-10 lived at `/tmp/workflow_drift_check.py` and reused `parse_groovydoc_full()` + `parse_includes()` from `bactopia-py/bactopia/nf.py`.

3. Present results as a clean summary:
   - **If all pass**: Report "All X modules/subworkflows have clean GroovyDoc" and stop
   - **If issues found**: Group by rule ID for easy batch fixing
     - Show count of affected components per rule
     - For each rule, list the affected components and specific messages
     - Prioritize FAIL over WARN

4. If the user asks to fix issues, apply fixes:

   **Module fixes:**
   - **M031 (FAIL)**: Read the actual `record()` output block, update `@output record(...)` to match
   - **M032 (FAIL)**: Read the actual input declaration, update `@input record(...)` to match
   - **M033 (FAIL)**: Reconcile optionality markers (`?`) between GroovyDoc and code. If the doc field has `?` but the code doesn't, either add `optional: true` / `Path?` to the code or drop the `?` from the doc. If the code is optional but the doc isn't, add the `?`.
   - **M034 (WARN)**: Remove description lines for standard fields (meta, results, logs, nf_logs, versions)
   - **M035 (FAIL)**: Check `data/citations.yml` for the correct key and update `@citation`
   - **M036 (WARN)**: Reorder tags to match: @status -> @keywords -> @tags -> @citation -> @note -> @input -> @output
   - **M037 (WARN)**: Add blank `*` lines between sections
   - **M038 (FAIL)**: Escape or reword any `*/` sequence inside the GroovyDoc block (e.g. replace glob patterns like `dir/*/file` with an alternative wording) so the comment doesn't close prematurely

   **Subworkflow fixes:**
   - **S017 (FAIL)**: Update `@modules` to match actual include statements (use underscore-delimited keys: `blast_blastn` not `blastn`)
   - **S018 (FAIL)**: Update `@subworkflows` to match actual include statements
   - **S019 (FAIL)**: Check `data/citations.yml` for the correct key and update `@citation`
   - **S023 (FAIL)**: Fix invalid feature values or add missing commas
   - **S024 (WARN)**: Reorder tags: @status -> @keywords -> @tags -> @citation -> @modules -> @subworkflows -> @note -> @input -> @output
   - **S026 (FAIL)**: Add `@output channel_name` for each undocumented emit channel, with field descriptions matching the emitted record shape (see pangenome subworkflow for pattern)
   - **S027 (FAIL)**: Move field descriptions from the `channel.empty()` @output to the correct @output that actually emits data

5. After fixes, re-run the lint to confirm all issues are resolved.

## Notes

- The wrapper script auto-discovers `bactopia-lint` (checks PATH, then conda envs)
- Citation keys are validated against `data/citations.yml` (all sections)
- Field ordering in @output is a WARN (style), not FAIL (correctness)
- Standard fields (meta, results, logs, nf_logs, versions) should never have description lines in @output
- The `@input` should use `record(meta, ...)` syntax, not raw `(_meta: Map, ...)` syntax

## GroovyDoc Quick Reference

### Required Tags

**Module order:** @status -> @keywords -> @tags -> @citation -> @note -> @input -> @output
**Subworkflow order:** @status -> @keywords -> @tags -> @citation -> @modules -> @subworkflows -> @note -> @input -> @output

```
@status stable|beta|deprecated
@keywords comma, separated, keywords
@tags complexity:<level> input-type:<type> output-type:<type> features:<list>
@citation comma, separated, bibtex_keys
```

### Subworkflow-Specific Tags
```
@modules bactopia_gather, csvtk_concat
@subworkflows scrubber, nohuman
```
Module names use underscore-delimited keys matching directory structure: `blast_blastn` (not `blastn`), `bactopia_qc` (not `qc`).

### Optional Tags
```
@note Title
Description text

@input record(meta, field1, field2)
- `meta`: Groovy Record containing sample information
- `field1`: Description

@input param_name
Description of non-record parameter

@output record(meta, field1, results, logs, nf_logs, versions)
- `field1`: Description (tool-specific fields only)
```

### Rules
- No spaces after commas in features list: `features:a,b,c` (not `features:a, b, c`)
- Standard fields (meta, results, logs, nf_logs, versions) must NOT have description lines
- @output field list must exactly match the actual `record()` output block
- @input record fields must match the actual input declaration
- @modules/@subworkflows must match actual include statements
- Citation keys must exist in `data/citations.yml`
- Valid complexity: simple, moderate, complex
- Valid input-type: none, single, multiple, parameter
- Valid output-type: single, multiple
