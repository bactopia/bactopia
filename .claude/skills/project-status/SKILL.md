---
name: project-status
description: Show a live snapshot of the Bactopia project state — component counts, GroovyDoc coverage, nf-test coverage, and structural issues. Use when asked about project state, coverage, what's missing, what's documented, or what needs attention.
---

# Project Status

Run the status script and interpret the output for the user.

## Steps

1. Run `bactopia-status` via the wrapper script:
   ```
   bash .claude/skills/project-status/scripts/run-bactopia-status.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --json
   ```

2. Parse the JSON output and present a clean, readable summary. Lead with:
   - Component counts (modules / subworkflows / workflows)
   - GroovyDoc coverage per tier (doc_count vs total)
   - Tag coverage per tier — highlight any required tags below 100%
   - nf-test coverage (workflows only)
   - Any structural issues (missing required files)

3. Highlight the most important gaps — prioritize by:
   - Missing required files (module.config, schema.json) first — these are structural breaks
   - Missing GroovyDoc second — documentation gaps
   - Missing nf-tests last — coverage gaps

4. If the user asks "what should I work on next?", suggest the highest-priority gaps from the output.

## Notes

- The wrapper script auto-discovers `bactopia-status` (checks PATH, then conda envs)
- **nf-test coverage is workflows-only** — modules and subworkflows are excluded for now
- Report whether the root pipeline test exists (`nftest.root_test`)

### JSON Output Fields
- `timestamp` — when the report was generated
- `repo` — path to the Bactopia repository
- `git` — `branch`, `commit`, `modified` file count
- `counts` — `modules`, `subworkflows`, `workflows`
- `tiers` — per-tier detail, each with `total`, `doc_count`, `issues[]`, `tag_coverage`
  - `tag_coverage` — map of tag name to count of components that have it
- `nftest` — `total` (workflow count), `tested`, `root_test`
### GroovyDoc Fields
All components use these required tags:
- `@status` — `stable`, `beta`, or `deprecated`
- `@keywords` — comma-separated search terms
- `@tags` — `complexity:<simple|moderate|complex> input-type:<single|multiple|parameter> output-type:<single|multiple> features:<list>`
- `@citation` — comma-separated bibtex keys
- `@note` — optional caveats or requirements

Modules and subworkflows also use:
- `@input` — input channel descriptions (can be multi-line for records)
- `@output` — output channel descriptions (single-line)

Subworkflows additionally use:
- `@subworkflows` — list of included subworkflows (directory keys)
- `@modules` — list of included modules (directory keys)

Entry workflows use:
- `@subworkflows` — included subworkflows (directory keys); workflows never call modules directly, so `@modules` is not used
- `@input` — parameter inputs (e.g., `rundir`)
- `@section` — groups published outputs
- `@publish` — file patterns with descriptions (within sections)
- `@note` — can appear at top level or within a section

### Required Files per Tier
- Modules: `module.config`, `schema.json`
- Subworkflows: none beyond `main.nf`
- Workflows: none beyond `main.nf`
