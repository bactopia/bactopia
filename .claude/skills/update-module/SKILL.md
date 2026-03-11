---
name: update-module
description: Check for newer versions of tools used in Bactopia modules and apply updates to module.config files and CHANGELOG.md. Use when asked to update module versions, check for outdated tools, or bump container versions.
---

# Update Module

Check for newer versions of bioconda tools used in Bactopia modules and apply updates.

## Steps

1. Run `bactopia-update` via the wrapper script:
   ```
   bash .claude/skills/update-module/scripts/run-bactopia-update.sh --bactopia-path /home/rpetit3/repos/bactopia/bactopia --json --silent
   ```
   If the user specified a module name, add `--module <name>` to the command.

2. Parse the JSON output. Separate entries into three categories:
   - **Needs update** (`needs_update: true`): ready for automatic update
   - **Needs user review** (`needs_user_review: true`): multi-package toolName, cannot auto-update
   - **Up to date** (`needs_update: false`): no action needed

3. Present results as a clean summary:
   - For updatable modules, show a table: tool, config path, installed version -> latest version
   - For user-review modules, list them separately with config path and a note that they require manual review
   - If nothing needs updating and nothing needs review, report "All modules are up to date" and stop

4. Ask the user which updates to apply (all or a subset). Do NOT proceed without confirmation.

5. For each confirmed update, read the module.config file at the path given by the `config` field, then edit three lines per module:

   **ext.toolName** -- Replace the conda spec string (the text between the first `"` and the `"` before `.replace`) with the `latest_toolName` value. Keep the `.replace()` chain untouched.
   ```
   // Before:
   ext.toolName = "bioconda::bakta=1.11.4".replace("=", "-").replace(":", "-").replace(" ", "-")
   // After (latest_toolName = "bioconda::bakta=1.12.0"):
   ext.toolName = "bioconda::bakta=1.12.0".replace("=", "-").replace(":", "-").replace(" ", "-")
   ```

   **ext.docker** -- Replace the entire quoted value with `latest_docker`.
   ```
   // Before:
   ext.docker = "biocontainers/bakta:1.11.4--pyhdfd78af_0"
   // After:
   ext.docker = "biocontainers/bakta:1.12.0--pyhdfd78af_0"
   ```

   **ext.image** -- Replace the entire quoted value with `latest_image`.
   ```
   // Before:
   ext.image = "https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0"
   // After:
   ext.image = "https://depot.galaxyproject.org/singularity/bakta:1.12.0--pyhdfd78af_0"
   ```

   For multi-process config files (e.g., genotyphi/parse has both genotyphi and mykrobe), use the `tool` field to match the correct `ext.toolName` line containing that tool name.

6. Update `CHANGELOG.md`. Under `## Unreleased Changes` > `### \`Added\``:
   - Look for an existing `- bump program versions in modules` line
   - If found, merge new entries into the indented list below it
   - If not found, add the block at the end of the `### \`Added\`` section
   - Format each entry as: `    - \`tool\`: old_version -> new_version`
   - One entry per tool (not per module) -- deduplicate shared tools (e.g., bakta_run and bakta_download both use bakta)
   - Sort entries alphabetically by tool name
   - Use the `tool` field (bioconda package name) as the display name

## Notes

- The wrapper script auto-discovers `bactopia-update` (checks PATH, then conda envs)
- Multi-package toolName modules (pirate, lissero, clonalframeml) have `needs_user_review: true` and no `latest_*` fields -- present them for manual review, do not auto-edit
- If `latest_version` is null for an entry (API failure), skip it and note the failure
- The `--module` flag filters by tool name prefix (e.g., `bakta` matches both bakta_run and bakta_download)

### JSON Output Fields

Standard update entry:
- `tool` -- bioconda package name
- `module` -- Bactopia module name (underscores)
- `config` -- relative path to module.config
- `installed_version` -- current version from module.config
- `latest_version` -- latest version from Anaconda API
- `latest_build` -- build string from Anaconda API
- `needs_update` -- true if versions differ
- `latest_toolName` -- new conda spec for ext.toolName
- `latest_docker` -- new docker image path for ext.docker
- `latest_image` -- new singularity image URL for ext.image

User review entry (multi-package):
- `tool`, `module`, `config`, `installed_version` -- same as above
- `needs_user_review` -- true (no latest_* fields provided)
