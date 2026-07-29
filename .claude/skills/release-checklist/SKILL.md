---
name: release-checklist
description: Audit whether Bactopia is ready for a version release and produce a GO / NO-GO recommendation report. Read-only — never commits, pushes, tags, edits issues, or regenerates tracked files. Use this whenever the user asks about release readiness, cutting a release, a release checklist, whether we can ship/tag a new version, pre-release audit, or "are we ready to release" — even if they don't say the word "checklist". Covers version consistency across bactopia/bactopia-py/nf-bactopia, nf-bactopia plugin pin + template drift, module tool-version bumps, workflow config/schema and catalog freshness, docs/citations/GroovyDoc/lint validators, recent test-run status, CHANGELOG completeness, and open-issue triage.
---

# Release Checklist

Audit whether Bactopia is ready to cut a version release and produce a recommendation report. This is a **multi-repo, audit-only** skill: it reads `bactopia`, `bactopia-py`, `nf-bactopia`, and the docs-site `bactopia.github.io` (published at https://bactopia.io), plus GitHub issues, and synthesizes a GO / NO-GO report.

## Hard constraints — read first

- **Audit-only. Never mutate anything.** No `git commit`/`push`/`tag`, no `gh issue edit`/`comment`/`close`, no edits to any tracked source file, no running the config/catalog/schema *regenerators* against the repo. The only write is the report file under `logs/` and scratch files under a `mktemp` dir.
- **Every regeneration check writes to a temp dir and diffs** — it never overwrites tracked `nextflow.config` / `nextflow_schema.json` / `catalog.json` / `llms.txt`.
- **Report, don't fix.** When a check finds drift, record it and recommend the follow-up skill (e.g. `/update-catalog`, `/merge-schemas`, `/update-module`). Do not apply the fix here.

## Version model (why the checks look where they do)

The Bactopia pipeline version and the `nf-bactopia` plugin pin are **declared in `versions.yml` at the pipeline repo root** (keys `bactopia:` and `nf-bactopia:`). `bactopia-merge-schemas` reads them and renders every workflow's `nextflow.config` from the bactopia-py template, and `bactopia-catalog` reads them for `catalog.json`, so the main repo's `nextflow.config` and `catalog.json` versions are **generated artifacts**. `versions.yml` is the version source of truth; bactopia-lint's repo rules (V001/V002/V003) enforce that every version-bearing file, the CHANGELOG heading, and the plugin pin agree with it. `CITATION.cff`, `bin/bactopia`, and `data/conda/meta.yaml` carry hand-maintained literals that must be bumped to match and are frequent laggards.

## Preamble

Resolve paths and a timestamp (from the bash tool):

```
BP=/home/rpetit3/repos/bactopia/bactopia
PY=$BP/../bactopia-py
NFB=$BP/../nf-bactopia
DOCS=$BP/../bactopia.github.io
TS=$(date +%Y%m%d_%H%M%S)
AUDIT=$(mktemp -d)          # scratch dir for read-only regen diffs
```

If a sibling repo path does not exist, run its checks in "skipped (repo not found)" mode and note it — never error out.

## Run the checks

Present each result; never fix. Numbers below are the report's check IDs.

### 1. Version consistency (delegated to bactopia-lint)

Version-bearing file consistency is enforced by bactopia-lint's **repo rules**, not re-implemented in this skill. Run bactopia-lint once — the same run feeds checks 8 and 11 — and read the `repo`-tier results:

```
bash .claude/skills/review-groovydoc/scripts/run-bactopia-lint.sh --bactopia-path $BP --json --silent
```

Parse the `repo` component's `results[]` for rule IDs starting with `V`:
- **V001** (FAIL): a version-bearing file disagrees with `versions.yml` — `nextflow.config`, `catalog.json`, `bin/bactopia`, `data/conda/meta.yaml`, `CITATION.cff`, or any `*.config` `bactopia_version` / `nf-bactopia@` pin. Module/subworkflow test configs inherit both from `conf/test_base.config`, so that single file is the usual test-side offender (not the ~180 individual test configs). The message groups offenders by value.
- **V002** (FAIL): `versions.yml: bactopia` ≠ the top `## vX.Y.Z` CHANGELOG heading (disagreement on what's being released).
- Passing rules are **not** listed individually (PASS components have empty `results[]`); absence of a V-FAIL means clean.
- Any V001/V002 FAIL ⇒ **check 1 = FAIL**. Fixes: regenerate configs (`/merge-schemas` + `/update-catalog`), hand-bump `CITATION.cff` / `bin/bactopia` / `data/conda/meta.yaml`, and bump `bactopia_version` in `conf/test_base.config` (the single source for every module/subworkflow test config).
- **Staleness guard:** the lint JSON MUST contain a `repo` component. If it does not, bactopia-py predates the version rules and the core version gate could not run — this is a **FAIL** (never a false PASS): report "version rules V001–V003 not loaded; update bactopia-py" and treat the release as NO-GO until re-audited with current tooling.

The deterministic engine is run once and feeds checks 2, 12, 13:

```
python3 .claude/skills/release-checklist/scripts/release_audit.py --bactopia-path $BP --json
```

It returns `sibling_release_state` (check 12), `changelog` (check 13), `module_updates` (check 2), and `warnings[]` — surface any warnings.

### 2. Module tool versions (offline — reads the /update-module record)

This check does **not** hit the network. `bactopia-update` queries the Anaconda API for ~100 modules (~2 min), so the audit never runs it — instead it reads the record that `/update-module` leaves under `logs/module-updates/<TS>.json`. The engine's `module_updates` field carries the result:

- `present: false` ⇒ **FAIL (blocking → NO-GO)**: no module-version check is on record for this cycle. Recommend: run `/update-module` (which writes the record) and then re-run `/run-tests`, before re-auditing. This is a hard release gate — do not hand-wave it.
- `present: true`:
    - `module_config_changed_after: true` ⇒ **WARN**: `module.config` files changed after the record was written, so it is stale — re-run `/update-module`.
    - `needs_update > 0` ⇒ **WARN**: that many modules have newer tool versions available — run `/update-module` to apply them.
    - otherwise ⇒ **PASS**: module versions were checked at `log_timestamp`; the `needs_user_review` count is informational (multi-package modules needing manual review).

Report `log_timestamp`, `needs_update`, `needs_user_review`, and `up_to_date` from the record. Never edit configs and never call `bactopia-update` from this skill.

### 3. Workflow config/schema freshness (read-only diff)

For each of the 4 named workflows — `bactopia`, `teton`, `staphopia`, `cleanyerreads` (paths from `catalog.json.workflows[<wf>].path`) — regenerate into the scratch dir and diff against the committed copies:

```
bash .claude/skills/merge-schemas/scripts/run-bactopia-merge-schemas.sh \
    --bactopia-path $BP --wf <wf> --outdir $AUDIT/<wf>
diff $AUDIT/<wf>/nextflow.config  <committed nextflow.config for wf>
diff $AUDIT/<wf>/nextflow_schema.json <committed nextflow_schema.json for wf>
```

The scratch dir is empty, so `--force` is not needed and nothing tracked is touched. Any diff ⇒ the committed file is stale relative to current module schemas/template ⇒ recommend `/merge-schemas`. The root `bactopia` config round-trips cleanly, so a root diff usually means a module was added/removed without rewiring, or the template version was bumped without regenerating. Only extend to the ~66 bactopia-tools if the user explicitly asks ("all tools") — tool schemas are otherwise covered by their nf-tests.

### 4. Catalog & llms.txt freshness (read-only diff)

```
bash .claude/skills/update-catalog/scripts/run-bactopia-catalog.sh \
    --bactopia-path $BP --output $AUDIT/catalog.json --pretty --llms-output $AUDIT/llms.txt
diff $AUDIT/catalog.json $BP/catalog.json
diff $AUDIT/llms.txt     $BP/llms.txt
```

Any diff ⇒ stale ⇒ recommend `/update-catalog`.

### 5. Docs-site (bactopia.io) version state

In `$DOCS` (skip with a note if the repo is absent):

- **(a) version label** — read `versions.current.label` in `docusaurus.config.ts`; compare to the release target. A lagging label is a WARN ("bump on release").
- **(b) changelog mirror** — compare the top `## v<x.y.z>` heading of `docs/changelog.md` to the main repo `CHANGELOG.md` top heading. A version mismatch is a WARN (mirror out of sync).
- **(c) plugin pins** — V001 only scans the bactopia repo's `*.config` files, so the docs site is not covered by the linter. Grep it directly: `grep -rn 'nf-bactopia@' "$DOCS" --include='*.md' --include='*.mdx'` and flag any pin that lags `versions.yml: nf-bactopia` (currently `developers/nf-bactopia/index.mdx` pins `2.0.3`). WARN.

### 6. Docs sync (.claude/docs)

```
bash .claude/skills/review-docs/scripts/run-bactopia-docs.sh --bactopia-path $BP --validate --json --silent
```

Record `summary.fail`. Any FAIL ⇒ recommend `/review-docs`.

### 7. Citations integrity

```
bash .claude/skills/review-citations/scripts/run-bactopia-citations.sh --bactopia-path $BP --validate --json --silent
```

Record `summary.orphans_total` and `summary.missing_total`. Either > 0 ⇒ recommend `/review-citations`. (`expected_orphans` are informational, not failures.)

### 8. GroovyDoc / lint

Reuse the bactopia-lint run from check 1 (do not run it twice). Count **component** FAIL results — module (`M0xx`/`MC0xx`), subworkflow (`S0xx`), workflow (`W0xx`) — and **exclude the `repo`-tier `V0xx` rules**, which are reported under checks 1 and 11. Any component FAIL ⇒ recommend `/review-groovydoc`.

### 9. Python lint (bactopia-py)

Run the configured ruff linter in `$PY` (prefer the justfile recipe, fall back to ruff directly):

```
cd $PY && just lint    # == poetry run ruff check .  ; fallback: ruff check .
```

Interpret the exit: lint findings (ruff reports violations) ⇒ **FAIL** recommending a ruff pass in bactopia-py. But distinguish "linter unavailable" from "lint failed": if neither `just` nor `ruff` is installed/resolvable in `$PY` (command-not-found, or `just lint` errors because the poetry env isn't set up), report **SKIP** ("ruff not available in this environment"), never FAIL — a missing tool is not a release blocker.

### 10. Test-run freshness (never launches a run)

Test runs live under `logs/run-tests/<timestamp>/` (written by `/run-tests`); `bactopia-review-tests` defaults to the newest one. An aborted run leaves a dir with **no `summary.json`** and the CLI errors on it, so pick the newest run dir that actually has a summary and review it explicitly:

```
LATEST_RUN=$(for d in $(ls -1dt "$BP"/logs/run-tests/[0-9]*/ 2>/dev/null); do [ -f "$d/summary.json" ] && basename "$d" && break; done)
bash .claude/skills/review-tests/scripts/run-bactopia-review-tests.sh --bactopia-path $BP --run "$LATEST_RUN" --silent
```

Also count how many `logs/run-tests/[0-9]*/` dirs are newer than `$LATEST_RUN` but lack a `summary.json` — those are incomplete/aborted runs worth flagging.

Relay the reviewed run's pass/fail status. Then judge **staleness**: compare `$LATEST_RUN`'s timestamp to the newest tracked-source commit (`git -C $BP log -1 --format=%ct`) and the working-tree dirty state (`sibling_release_state.bactopia.dirty_files` from the engine). Report `stale` if code changed after the run or the tree is dirty. Severity: a failing reviewed run ⇒ **FAIL**; a passing-but-stale run, or the presence of newer incomplete run dirs, ⇒ **WARN** recommending `/run-tests`; if no run dir has a `summary.json` at all ⇒ **WARN** ("no completed test run on record"). Never trigger a run from this skill.

### 11. nf-bactopia plugin currency (bactopia-lint V003)

From the same lint run, read repo rule **V003**: `versions.yml: nf-bactopia` vs the nf-bactopia repo's `build.gradle` latest. A `WARN` means the declared pin lags the newest nf-bactopia release — surface it; adopting the newer plugin may be intentional, so it's a judgment call, not a hard block. PASS/absent = current. Per-config pin drift is already covered by V001 (check 1). If the nf-bactopia repo isn't checked out, V003 PASSes with a "cannot verify" note.

### 12. Sibling repo release state

From the engine's `sibling_release_state`: per repo (`bactopia`, `bactopia-py`, `nf-bactopia`) report `latest_tag`, `commits_ahead`, `dirty_files`, `changelog_top`, `needs_release`. Flag each repo with `needs_release: true` (WARN) — unreleased commits past its last tag mean it may need its own release before the pipeline release.

### 13. CHANGELOG completeness & concision

Model-judgment check using the engine's `changelog` data:

- `placeholder: true` (codename `"???"` or a `?` in the date) ⇒ **FAIL** (fill codename + date before release).
- `missing_section_for_target: true` ⇒ **FAIL** (the declared `target_version` from `versions.yml` has no `## v<target>` CHANGELOG section at all). Whether the top heading equals the target is V002's job (check 1).
- Concision: the last git tag can predate multiple unreleased versions, so `commits_since_tag` may legitimately exceed `top_section_bullets` — do **not** treat that as an equality target. WARN only if top-section bullets are multi-sentence/verbose relative to peers, or the coverage gap looks like genuinely missing entries. Keep it qualitative; no bullet cap.

### 14. Open-issue triage

```
gh issue list --repo bactopia/bactopia --state open --limit 300 \
    --json number,title,labels,updatedAt
```

- Exclude issues whose labels include the ignore label — default `release-ignore`, or whatever label the user names. If the label doesn't exist in the repo, nothing is filtered (no error). Report how many were excluded.
- Report total open, the excluded count, and a breakdown by label.
- Surface a **shortlist (≤10)** of likely release-relevant issues: labeled `bug`, referenced by `#<number>` in the top CHANGELOG section, or updated within ~30 days. For each, give a one-line **address / postpone** suggestion.
- This check is **INFO** — it never changes the overall verdict on its own. Never edit, label, or comment on issues.

### 15. Version-pinned datasets published

`bactopia datasets` downloads pre-compiled bundles from `https://datasets.bactopia.com/datasets/v<version>/`, where `<version>` is the pipeline version from `versions.yml`. On a release the version is bumped **before** that version's bundle is uploaded, so the URL 404s and every dataset-dependent run — and the `bactopia_datasets` nf-test — fails. Probe the target version's bundle (read-only HEAD; the audit never uploads):

```
V=$(awk '/^bactopia:/{print $2}' $BP/versions.yml)
curl -s -o /dev/null -w '%{http_code}' -I "https://datasets.bactopia.com/datasets/v${V}/amrfinderplus.tar.gz"
```

- HTTP `200` ⇒ **PASS**: the v`<version>` datasets are published.
- Anything else (typically `404`) ⇒ **FAIL (blocking → NO-GO)**: the v`<version>` dataset bundle is not published. Publishing the versioned datasets is a required release step; until it lands, `bactopia datasets` 404s for users and the `bactopia_datasets` tests fail. The version-pinned `amrfinderplus.tar.gz` is a sufficient sentinel (other bundles like `mash-refseq88...` are not version-gated). **Fix:** run `/update-datasets` to rebuild `amrfinderplus.tar.gz` in the pinned container and publish it to `bactopia-r2:bactopia/datasets/v<version>/`.

## Synthesize the report (the primary deliverable)

The report is produced two ways from the same content:

1. **Write the full report** to `$BP/logs/release-audit-$TS.md` using the template below (all 15 checks, Blocking / Non-blocking sections, and a Details section per check with findings + the recommended follow-up skill). `logs/` is the scratch/output dir; use a single timestamped **file** (not a `logs/<TS>/` dir) so `/review-tests`' run-dir scan ignores it. If `logs/` is unwritable, fall back to `$AUDIT/release-audit-$TS.md` and report that path.
2. **Present a condensed view in chat**: the Overall verdict line, the Checklist table, and the Blocking-items list — then end with the saved path, e.g. `Full report saved to: /home/rpetit3/repos/bactopia/bactopia/logs/release-audit-<TS>.md`.

### Report template

```
# Bactopia Release Readiness Audit — <TS>
**Target release:** v<versions.yml: bactopia>   **Overall: <GO | GO WITH CAVEATS | NO-GO>**

<one-paragraph rationale>

## Checklist
| # | Check | Status | Summary |
|---|-------|--------|---------|
| 1  | Version consistency          | <PASS/WARN/FAIL> | ... |
| 2  | Module tool versions         | ... | ... |
| 3  | Workflow configs & schemas   | ... | ... |
| 4  | Catalog & llms.txt freshness | ... | ... |
| 5  | Docs-site (bactopia.io) state| ... | ... |
| 6  | Docs sync (.claude/docs)     | ... | ... |
| 7  | Citations                    | ... | ... |
| 8  | GroovyDoc / lint             | ... | ... |
| 9  | Python lint (ruff)           | ... | ... |
| 10 | Test-run freshness           | ... | ... |
| 11 | nf-bactopia plugin currency  | ... | ... |
| 12 | Sibling repo release state   | ... | ... |
| 13 | CHANGELOG completeness       | ... | ... |
| 14 | Open issues                  | INFO | <N open, M excluded> |
| 15 | Version-pinned datasets      | <PASS/FAIL> | <v{version} bundle: 200 / 404> |

## Blocking items (must fix before release)
- ...

## Non-blocking / judgment items
- ...

## Details
### 1. Version consistency (bactopia-lint V001/V002)
<V001/V002 findings; fix = run `/bump-versions` (propagates versions.yml -> conf/test_base.config + CITATION.cff/bin/bactopia/meta.yaml), then regenerate generated artifacts with /merge-schemas + /update-catalog>
### 2. Module tool versions
...
```

### Severity + overall verdict (deterministic — status is never a guess)

- **FAIL (blocking)**: check 1 — bactopia-lint **V001** (a version-bearing file ≠ `versions.yml`) or **V002** (`versions.yml` ≠ CHANGELOG top heading) reports FAIL, **or** the lint `repo` component is absent (version rules not loaded); CHANGELOG `placeholder` or `missing_section_for_target` (check 13); **no `/update-module` record — `module_updates.present == false` (check 2)**; schema/config drift (check 3); catalog/llms drift (check 4); any component docs/citations/lint failure (checks 6–8); a ruff lint failure when ruff is available (check 9); a failing test run (check 10); **the target-version dataset bundle is unpublished — `datasets.bactopia.com/datasets/v<version>/` 404s (check 15)**.
- **WARN (non-blocking / judgment)**: pending module updates or a stale `/update-module` record (check 2); the declared nf-bactopia pin lags latest — V003 (check 11); a sibling repo with `needs_release` (check 12); docs-site label or changelog-mirror lag (checks 5a/5b); a passing-but-stale test run (check 10); a CHANGELOG concision concern (check 13).
- **PASS / SKIP / INFO**: PASS = check clean; SKIP = a check that could not run (e.g. ruff not installed, check 9) — never a blocker; INFO = informational only (check 14). None of these force a caveat.
- **Overall** = `NO-GO` if any FAIL; `GO WITH CAVEATS` if only WARN; `GO` if all PASS. Open issues (check 14) are INFO and never change the overall verdict on their own.

## Notes

- All `run-bactopia-*.sh` wrappers auto-discover their CLI (PATH → `bactopia-dev` conda env → `bactopia-py` → any `bactopia-*` env), so no env activation is needed.
- `--bactopia-path` is always `/home/rpetit3/repos/bactopia/bactopia`.
- The engine (`release_audit.py`) needs only `git` + the Python stdlib; it does not use conda or any `bactopia-*` CLI.
