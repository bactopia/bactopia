# Skills Reference

## Overview

Skills are AI tooling — short instruction files that Claude invokes via the `Skill` tool when a matching trigger phrase appears. They are distinct from pipeline components (modules, subworkflows, workflows) and are **not** listed in [catalog.json](../../../catalog.json); this doc is the authoritative inventory for project-local skills.

**Pattern**: every project-local skill is a thin AI wrapper over a `bactopia-*` CLI in [bactopia-py](../../../../bactopia-py/bactopia/cli/). The skill lives at `.claude/skills/<name>/SKILL.md`; a wrapper script under `scripts/` discovers and invokes the CLI. Implementation logic belongs in the CLI so it can evolve independently and be run directly from the shell or CI — the `SKILL.md` is just interpretation.

## Project-local skills

| Skill | Backend | Purpose |
|---|---|---|
| [add-bactopia-tool](../../skills/add-bactopia-tool/) | `bactopia-scaffold` | Scaffold a complete Bactopia Tool across all three tiers -- module, subworkflow, and workflow entry point under workflows/bactopia-tools/. |
| [add-module](../../skills/add-module/) | `bactopia-scaffold` | Scaffold a new Bactopia module from a bioconda/conda-forge package. |
| [add-subworkflow](../../skills/add-subworkflow/) | `bactopia-scaffold` | Scaffold a new Bactopia subworkflow that orchestrates existing modules. |
| [merge-schemas](../../skills/merge-schemas/) | `bactopia-merge-schemas` | Regenerate nextflow.config and nextflow_schema.json for Bactopia workflows by running bactopia-merge-schemas. |
| [project-status](../../skills/project-status/) | `bactopia-status` | Show a live snapshot of the Bactopia project state — component counts, GroovyDoc coverage, nf-test coverage, and structural issues. |
| [review-citations](../../skills/review-citations/) | `bactopia-citations` | Review citation integrity across data/citations.yml and @citation tags using bactopia-citations --validate. |
| [review-docs](../../skills/review-docs/) | `bactopia-docs` | Review staleness of reference docs under .claude/docs/ using bactopia-docs --validate. |
| [review-groovydoc](../../skills/review-groovydoc/) | `bactopia-lint` | Review GroovyDoc accuracy across modules and subworkflows using bactopia-lint. |
| [review-tests](../../skills/review-tests/) | — | Review nf-test run results and present a diagnostic summary with grouped error analysis. |
| [run-tests](../../skills/run-tests/) | `bactopia-test` | Run Bactopia nf-tests via bactopia-test and produce a timestamped logs/ directory that /review-tests can interpret. |
| [update-catalog](../../skills/update-catalog/) | `bactopia-catalog` | Regenerate catalog.json and llms.txt by running bactopia-catalog. |
| [update-module](../../skills/update-module/) | `bactopia-update` | Check for newer versions of tools used in Bactopia modules and apply updates to module.config files and CHANGELOG.md. |

The `Purpose` column is the first sentence of each skill's `description:` frontmatter. Drift between the table and the source file is caught by **D107** in `/review-docs`. Full trigger-phrase lists live in each `SKILL.md` — read it directly when you need the exact phrasing.

## Global skills worth knowing

These live at `~/.claude/skills/` (user-global) or are built into the Claude Code harness:

| Skill | When to use |
|---|---|
| `skill-creator` | **Any time** you are creating, editing, or improving a skill. Do not hand-scaffold `SKILL.md` files. |
| `session-handoff` | Save session context to the clipboard for resuming in a new session. |
| `update-config` | Configure the Claude Code harness via `settings.json` (hooks, permissions, status line). |
| `schedule` / `loop` | Recurring cron-style agents (`schedule`) or self-paced polling (`loop`). |
| `claude-api` | Reference for building apps on the Anthropic SDK — unrelated to the Bactopia pipeline. |

## When to invoke `skill-creator`

Use it whenever the request is about *making* or *modifying* a skill:

- "Create a new skill that does X" / "Add a `/foo` slash command"
- "Edit the `<name>` skill to also do Y"
- "Improve the `<name>` skill" / "Make `<name>` handle Z"
- "Scaffold a skill for …"

`skill-creator` owns the conventions for frontmatter, naming, trigger-phrase design, and output formatting. Skills created outside it will usually need to be reworked to match — use it up front.

## Architecture

```text
.claude/skills/<name>/
    SKILL.md                # Frontmatter (name, description) + trigger instructions
    scripts/                # Optional — wrapper scripts
        run-bactopia-*.sh   # Discovers the bactopia-* CLI in PATH or conda envs
```

Wrapper scripts find the backing CLI in this order: PATH first, then the exact bactopia-dev conda env, then exact bactopia-py, then any `bactopia-*` env (fuzzy match). This discovery code is shared across every wrapper, so env-setup changes propagate uniformly.

**Adding a new CLI-backed skill**: implement the CLI in [bactopia-py](../../../../bactopia-py/bactopia/cli/) first so it's testable, stable, and runnable outside Claude, then wrap it with `skill-creator`. The skill becomes the trigger surface; the CLI is the logic.

## Discoverability

- **Terminal CLI**: custom skills appear in the `/` autocomplete menu.
- **VS Code extension**: custom skills do **not** appear in the slash-command menu (known limitation). They still work — type `/skill-name` directly or describe what you want and Claude will trigger them via the same system reminder that enumerates skills at session start.
- This reference doc is linked from [CLAUDE.md](../../../CLAUDE.md) and [llms.txt](../../../llms.txt) so Claude has a canonical inventory to consult.
