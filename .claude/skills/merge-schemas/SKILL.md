---
name: merge-schemas
description: Regenerate nextflow.config and nextflow_schema.json for Bactopia workflows by running bactopia-merge-schemas. Use when asked to merge schemas, regenerate workflow config, rebuild nextflow_schema.json, or sync workflow configs after module schema changes.
---

# Merge Schemas

Regenerate the `nextflow.config` and `nextflow_schema.json` files for one or more Bactopia workflows by invoking `bactopia-merge-schemas`. Auto-discovers each workflow's output directory from `data/catalog.json` so the caller never hand-computes paths.

## Steps

1. **Identify the target workflow(s).** Parse the user's request into one of these shapes:
    - A specific workflow name (e.g. `teton`, `merlin`, `bactopia`)
    - A list (e.g. `teton and staphopia`) or a glob-ish pattern (e.g. `all merlin tools`, `everything starting with blast`)
    - `all tools` → every workflow with `type: tool` in the catalog
    - `all named` → every workflow with `type: named`
    - `all` → every workflow in the catalog

    If the ask is ambiguous (e.g. "regenerate merlin" — merlin is a named Bactopia Tool entry but also a subworkflow tag), confirm with the user before running anything. Never run without some kind of filter.

2. **Resolve each target to a path from the catalog.** Use the Read tool to read `/home/rpetit3/repos/bactopia/bactopia/data/catalog.json`. For each requested workflow `W`:
    - Verify `catalog["workflows"][W]` exists. If not, stop and report the unknown name — do not invoke the CLI.
    - Pull `catalog["workflows"][W]["path"]` — this is the exact value to pass as `--outdir`. It is already relative to the repo root:
        - `.` → main `bactopia` workflow (outputs land in the repo root)
        - `workflows/teton/` / `workflows/staphopia/` / `workflows/cleanyerreads/` → named workflows
        - `workflows/bactopia-tools/{name}/` → Bactopia Tools
    - Record `catalog["workflows"][W]["type"]` (`named` vs `tool`) for the summary.

3. **Check for existing generated files.** For each target `--outdir`, use Read or Glob to check whether `nextflow_schema.json` or `nextflow.config` already exist. If any do:
    - List the affected workflows and files for the user
    - Get explicit confirmation before proceeding with `--force`
    - Without `--force`, the CLI exits with an error on the first existing file

    Do not pass `--force` preemptively — it is a destructive overwrite of tracked files.

4. **Run `bactopia-merge-schemas` once per target** via the wrapper:
    ```
    bash .claude/skills/merge-schemas/scripts/run-bactopia-merge-schemas.sh \
        --bactopia-path /home/rpetit3/repos/bactopia/bactopia \
        --wf <workflow> \
        --outdir <path-from-catalog> \
        --force    # only if the user confirmed in step 3
    ```
    Run sequentially. On the first failure, stop and report — do not silently continue through a batch.

5. **Summarize the result.** For each target, report:
    - Workflow name and type (`named` / `tool`)
    - The `--outdir` used
    - Which files were generated (`nextflow_schema.json`, `nextflow.config`)
    - Pass/fail status

    If nothing was regenerated (e.g. all targets refused `--force`), say so explicitly.

## Notes

- The wrapper auto-discovers `bactopia-merge-schemas` (PATH → `bactopia-dev` conda env → `bactopia-py` conda env → any `bactopia-*` env)
- `--bactopia-path` is always `/home/rpetit3/repos/bactopia/bactopia` — do not guess or prompt
- `--outdir` always comes from `catalog["workflows"][wf]["path"]` — never hand-compute it
- The CLI writes only **two** files: `nextflow_schema.json` and `nextflow.config`. The source also declares `params.config` / `process.config` variables and includes them in its pre-flight existence check, but does not actually write them. If the user expects those files, tell them the current CLI does not generate them
- `--force` is required if *either* of the two output files already exists. Ask first
- `conf/schema/{wf}.json` backs named workflows; `conf/schema/bactopia-tools.json` backs tools. Generic parameters come from `conf/schema/generic.json`. You do not need to touch these files — the CLI reads them directly
- Typical catalog contents: 4 named workflows (`bactopia`, `teton`, `staphopia`, `cleanyerreads`) and ~66 tools under `workflows/bactopia-tools/`

### CLI Reference (`bactopia-merge-schemas`)

Required:
- `--bactopia-path PATH` — repo root
- `--wf NAME` — workflow key from `data/catalog.json`

Output:
- `--outdir PATH` — directory to write outputs (default `.`)
- `--force` — overwrite existing output files

Other:
- `--verbose` / `--silent` — log level
- `--version` / `--help`

### Sibling Skills

- `/update-module` — bumps tool versions in `module.config` and `CHANGELOG.md`. When a module update changes that module's `schema.json`, `/merge-schemas` is the right follow-up for every workflow that includes the module.
- `/project-status` — shows coverage and structural issues; useful for confirming catalog state before batch regeneration.
