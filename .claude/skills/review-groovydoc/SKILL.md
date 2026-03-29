---
name: review-groovydoc
description: Review GroovyDoc accuracy across all modules using bactopia-lint. Checks @output/@input field matching, citation keys, tag ordering, and formatting. Use when asked to review GroovyDoc, check documentation accuracy, validate module docs, or audit GroovyDoc.
---

# Review GroovyDoc

Run bactopia-lint focused on GroovyDoc accuracy rules and present the results.

## Steps

1. Run `bactopia-lint` via the wrapper script with `--modules --quiet --json --silent`:
   ```
   bash .claude/skills/review-groovydoc/scripts/run-bactopia-lint.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --modules --no-subworkflows --no-workflows --quiet --json --silent
   ```
   If the user specified a module name, add `--module <name>` to the command.

2. Parse the JSON output. Filter results to **GroovyDoc-related rules only**:
   - M006: GroovyDoc block present
   - M007: Required GroovyDoc tags
   - M008: @status value valid
   - M009: @tags sub-keys present
   - M010: Features list formatting
   - M014: Links use HTTPS
   - M031: @output fields match actual record() output
   - M032: @input record fields match actual input
   - M033: @input uses record() syntax
   - M034: Standard fields not described in @output
   - M035: @citation keys exist in citations.yml
   - M036: GroovyDoc tag ordering
   - M037: Blank lines between GroovyDoc sections

3. Present results as a clean summary:
   - **If all pass**: Report "All X modules have clean GroovyDoc" and stop
   - **If issues found**: Group by rule ID for easy batch fixing
     - Show count of affected modules per rule
     - For each rule, list the affected modules and specific messages
     - Prioritize FAIL over WARN

4. If the user asks to fix issues, apply fixes:
   - **M031 (FAIL)**: Read the actual `record()` output block, update `@output record(...)` to match
   - **M032 (FAIL)**: Read the actual input declaration, update `@input record(...)` to match
   - **M033 (WARN)**: Convert `@input (_meta: Map, ...)` to `@input record(meta, ...)`
   - **M034 (WARN)**: Remove description lines for standard fields (meta, results, logs, nf_logs, versions)
   - **M035 (FAIL)**: Check `data/citations.yml` for the correct key and update `@citation`
   - **M036 (WARN)**: Reorder tags to match: @status -> @keywords -> @tags -> @citation -> @note -> @input -> @output
   - **M037 (WARN)**: Add blank `*` lines between sections

5. After fixes, re-run the lint to confirm all issues are resolved.

## Notes

- The wrapper script auto-discovers `bactopia-lint` (checks PATH, then conda envs)
- Citation keys are validated against `data/citations.yml` (all sections)
- Field ordering in @output is a WARN (style), not FAIL (correctness)
- Standard fields (meta, results, logs, nf_logs, versions) should never have description lines in @output
- The `@input` should use `record(meta, ...)` syntax, not raw `(_meta: Map, ...)` syntax

## GroovyDoc Quick Reference

### Required Tags (in order)
```
@status stable|beta|deprecated
@keywords comma, separated, keywords
@tags complexity:<level> input-type:<type> output-type:<type> features:<list>
@citation comma, separated, bibtex_keys
```

### Optional Tags
```
@note Title
Description text

@input record(meta, field1, field2)
- `meta`: Groovy Map containing sample information
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
- Citation keys must exist in `data/citations.yml`
