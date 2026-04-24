---
name: update-catalog
description: Regenerate catalog.json and llms.txt by running bactopia-catalog. Use when asked to update the catalog, rebuild the component index, refresh catalog.json, or sync llms.txt after component changes (new/removed modules, subworkflows, or workflows; tool version bumps; GroovyDoc edits that affect descriptions or contracts).
---

# Update Catalog

Regenerate the machine-readable Bactopia component index (`catalog.json`) and the AI-discovery surface (`llms.txt`) by invoking `bactopia-catalog`. Both files live at the repo root. The skill is a thin wrapper — all scanning, parsing, and rendering logic lives in `bactopia-catalog` in bactopia-py.

## Steps

1. **Tell the user** you are checking for uncommitted local edits to `catalog.json` and `llms.txt`, then run:
    ```
    git -C /home/rpetit3/repos/bactopia/bactopia status --porcelain catalog.json llms.txt
    ```
    If either file is dirty, show the user the diff and confirm before overwriting. `bactopia-catalog` has no `--force` flag and silently clobbers existing output files, so this check is the user's only safety net against losing in-flight edits. If neither file is dirty, proceed directly to Step 2.

2. **Run the wrapper** to regenerate both files in one invocation:
    ```
    bash .claude/skills/update-catalog/scripts/run-bactopia-catalog.sh \
        --bactopia-path /home/rpetit3/repos/bactopia/bactopia \
        --output /home/rpetit3/repos/bactopia/bactopia/catalog.json \
        --pretty \
        --llms-output /home/rpetit3/repos/bactopia/bactopia/llms.txt
    ```
    - Always pass `--pretty` so the committed `catalog.json` stays diff-friendly.
    - `--llms-output` is what triggers `llms.txt` regeneration. Without it, only `catalog.json` is written. The user asked to update both — always pass it unless they explicitly say "catalog only".
    - The default llms.txt template is bundled inside bactopia-py at `bactopia/templates/bactopia/llms.txt.j2`. Do not pass `--llms-template` unless the user asks you to render from a non-default template path.

3. **Summarize the result.** Use the CLI's own summary line for counts (`Found N modules, M subworkflows, K workflows`) — don't parse the JSON yourself. Then show the user:
    - `catalog.json` path and the shape of the change (`git diff --shortstat catalog.json`)
    - `llms.txt` path and whether anything actually changed (`git diff --shortstat llms.txt`)
    - If a count changed (new or removed module/subworkflow/workflow), suggest `/project-status` for a full coverage read.

4. **Do not stage or commit** the updated files. Let the user review the diff first. Committing auto-generated files without human review is how stale or malformed content sneaks into main.

## Notes

- **Wrapper discovery order**: PATH → `bactopia-dev` conda env → `bactopia-py` conda env → any `bactopia-*` env. Matches the shared pattern used by `merge-schemas`, `project-status`, `update-module`, and `run-tests`.
- **`--bactopia-path` is always `/home/rpetit3/repos/bactopia/bactopia`** — do not guess or prompt. This is a single-repo project.
- **File locations**: `catalog.json` and `llms.txt` both live at the repo root as of the catalog-location migration. Older references to `data/catalog.json` are stale.
- **`llms.txt` is auto-generated** — prose lives in the Jinja2 template at `bactopia-py/bactopia/templates/bactopia/llms.txt.j2` (inside the bactopia-py package, not in this repo). If the user wants to change prose (headings, category highlights, key-patterns bullets), they should edit the template in bactopia-py, not `llms.txt` here. Running this skill without first editing the template will regenerate a byte-identical `llms.txt`.
- **What's dynamic in the template**: the module count line, the "Workflows (Tier 1)" list (looped over named workflows from the catalog), and the tool count. Everything else is hand-curated prose in the template.
- **No `--force` guard in the CLI**: `bactopia-catalog` silently overwrites whatever path `--output` / `--llms-output` point at. Step 1's dirty-check is the pragmatic substitute — don't skip it.
- **Right follow-up after**: `/add-module`, `/add-subworkflow`, `/update-module` (when schema or tool version changes affect catalog contents), or any manual edit to a module/subworkflow/workflow `main.nf` that changes GroovyDoc descriptions, contracts, or tags.

### CLI Reference (`bactopia-catalog`)

Required:
- `--bactopia-path PATH` — directory where the Bactopia repository is stored

Catalog output:
- `-o, --output PATH` — where to write `catalog.json` (default: stdout, which is almost never what you want from a skill)
- `--pretty` — pretty-print JSON with 2-space indentation (always pass this so the committed file diffs cleanly)

llms.txt output:
- `--llms-output PATH` — where to write the rendered `llms.txt`. Omit to skip llms.txt rendering.
- `--llms-template PATH` — Jinja2 template path. Defaults to the bundled template inside bactopia-py (`bactopia/templates/bactopia/llms.txt.j2`). Rarely needed.

Other:
- `--verbose`, `--version`, `--help`

### Sibling Skills

- `/project-status` — reads `catalog.json` via `bactopia-status`. A freshly regenerated catalog is immediately reflected in its output, so this is the natural read-after-write check.
- `/merge-schemas` — also reads `catalog.json` (for workflow path resolution). Run this after `/update-catalog` when a workflow's module set or schema has changed.
- `/add-module`, `/add-subworkflow`, `/update-module` — all three mutate the component graph. `/update-catalog` is the standard follow-up to keep `catalog.json` and `llms.txt` in sync.
