---
name: bump-versions
description: Propagate the Bactopia and nf-bactopia versions declared in versions.yml into the hand-maintained files that carry a literal version (conf/test_base.config, CITATION.cff, bin/bactopia, data/conda/meta.yaml). Use this whenever the user has edited versions.yml and wants the rest of the repo brought in line, or asks to bump the version, set the release version, propagate versions.yml, sync version-bearing files, fix a V001 version-consistency failure, or prepare version files for a release. This never edits versions.yml itself and never regenerates the templated artifacts (nextflow.config, catalog.json) — it hands those off to /merge-schemas and /update-catalog.
---

# Bump Versions

`versions.yml` at the repo root is the **single source of truth** for the pipeline
version (`bactopia`) and the plugin pin (`nf-bactopia`). It is edited **by hand** —
this skill never touches it. Its job is the mechanical, error-prone part: copying
those two values into every *hand-maintained* file that repeats them, so
bactopia-lint's **V001** (version-bearing files must equal `versions.yml`) goes green.

Two classes of version-bearing files exist, and only one is this skill's concern:

- **Hand-maintained (this skill fixes these):** `conf/test_base.config`
  (`bactopia_version` + the `nf-bactopia@` plugin pin, inherited by all ~183
  module/subworkflow test configs), `CITATION.cff`, `bin/bactopia`,
  `data/conda/meta.yaml`.
- **Generated (this skill does NOT touch these):** `nextflow.config`, every
  `workflows/*/nextflow.config`, `catalog.json`, `llms.txt`. These are rendered
  from `versions.yml` by `bactopia-merge-schemas` / `bactopia-catalog`, so the fix
  is to *regenerate* them — a follow-up, not an edit. Slash commands cannot invoke
  other slash commands, so you recommend those to the user rather than running them.

## Steps

1. **Dry-run first.** Show the user exactly what will change before writing anything:
   ```
   python3 .agents/skills/bump-versions/scripts/bump_versions.py \
       --bactopia-path /home/rpetit3/repos/bactopia/bactopia --check
   ```
   The script reads `versions.yml`, then reports each hand-maintained literal as
   `would update <old> -> <new>`, `already <value>`, or a `WARN` (file missing or
   the literal moved). If everything is already `ok`, tell the user the
   hand-maintained files are in sync and skip to step 4 (they may still need a
   regen if a generated artifact lagged).

2. **Confirm, then apply.** Once the user is happy with the dry-run, drop `--check`
   to write the changes:
   ```
   python3 .agents/skills/bump-versions/scripts/bump_versions.py \
       --bactopia-path /home/rpetit3/repos/bactopia/bactopia
   ```
   Only the version token inside each match is rewritten; runs are idempotent.

3. **Recommend the regeneration follow-ups.** If anything changed (or a generated
   artifact is suspected stale), tell the user to run, in order:
   - `/merge-schemas all` — re-renders `nextflow.config` + workflow configs + schemas from `versions.yml`.
   - `/update-catalog` — rebuilds `catalog.json` + `llms.txt`.
   Do not attempt to run these from here.

4. **Verify V001–V003.** Confirm the version gate is actually green:
   ```
   bash .agents/skills/review-groovydoc/scripts/run-bactopia-lint.sh \
       --bactopia-path /home/rpetit3/repos/bactopia/bactopia --json --silent
   ```
   Read the `repo` component's results. A clean run lists **no** `V0xx` failures
   (PASS results are suppressed). If **V001** still fails, report the offending
   files verbatim — most often a generated artifact whose regen (step 3) hasn't run
   yet, or a version-bearing file this skill doesn't own.

5. **Flag downstream checks (report, don't fix):**
   - **V002** compares `versions.yml` to the top `## vX.Y.Z` CHANGELOG heading. If
     they disagree, the CHANGELOG needs a section for the new version — that's
     human-authored content, so surface it, don't invent it.
   - Bumping `bactopia_version` changes the test container tag
     (`bactopia/bactopia:<version>`) used by every module/subworkflow test. Recommend
     a smoke test (`/run-tests` on a quick component) once the matching image exists,
     since a nonexistent tag will fail tests even with a green lint.

## Notes

- The script needs only Python stdlib + read access to `versions.yml`; no conda env
  or `bactopia-*` CLI. There is deliberately **no** `bactopia-versions` CLI —
  `versions.yml` is hand-edited and this propagation is the whole story.
- The four hand-maintained targets and their literal formats are encoded in the
  script's `TARGETS` map, matching V001's own regexes so the skill fixes exactly what
  the linter checks. If V001 ever grows a new hand-maintained file, add it there.

### Sibling skills

- `/merge-schemas`, `/update-catalog` — the regeneration follow-ups (step 3).
- `/review-groovydoc` — backs the `bactopia-lint` verification run (step 4).
- `/release-checklist` — its check 1 delegates the version gate to these same V-rules.
