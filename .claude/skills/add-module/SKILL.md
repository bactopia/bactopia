---
name: add-module
description: Scaffold a new Bactopia module from a bioconda/conda-forge package. Creates main.nf, module.config, schema.json, and test files following project standards. Use when asked to add a new module, create a new module, scaffold module files, or add a new tool's module.
---

# Add Module

Scaffold a complete Bactopia module for a bioconda/conda-forge package, creating all required files with GroovyDoc documentation and nf-test tests.

## Prerequisites

Before using this skill, read:
- `.claude/docs/standards/05-module-documentation.md` -- Module standards including module.config, schema.json, and test templates
- `.claude/docs/project/04-testing-framework.md` -- Testing framework details

## Phased Workflow

Follow these phases in order. When unsure about ANYTHING, ask the user rather than guess.

---

### Phase 1: Package Verification

**Goal:** Confirm the bioconda package exists and retrieve version/container information.

1. Ask the user for the **bioconda package name** (e.g., `mlst`, `bakta`, `ssuissero`).

2. Ask: **Is this a standalone module or part of a multi-process set?**
   - Standalone: files go in `modules/{tool}/`
   - Multi-process: files go in `modules/{tool}/{process}/` (e.g., `modules/bakta/run/`)
   - If multi-process: ask which process this is (run, download, summary, collate, etc.)

3. Run the lookup command:
   ```bash
   bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-scaffold.sh lookup {package_name} --bactopia-path . --pretty
   ```

4. The output includes:
   - `package`, `channel`, `version`, `build` -- package identity
   - `summary`, `home` -- tool description and documentation URL
   - `container_refs` -- `toolName`, `docker`, `image` strings
   - `existing_components.module` -- whether the module already exists

5. **Present findings to the user** and ask them to confirm before proceeding.
   If the module already exists, warn the user.

---

### Phase 2: Information Gathering

**Goal:** Determine tool behavior, inputs, outputs, and parameters through web lookup + user interview.

1. **Fetch the tool's documentation** using WebFetch on the `home` URL from Phase 1.
   - If the URL is a GitHub repo, fetch the main page (includes README)
   - Extract: command-line options, input file types, output files, version command
   - If WebFetch fails, skip to step 2 and ask the user directly

2. **Ask the user to confirm or modify** the following (present what you found from docs):

   **a. Input type** -- "What does this tool take as input?"
   | User says | Input record fields |
   |-----------|---------------------|
   | Assembly / contigs / FASTA | `record(meta: Record, fna: Path)` |
   | Reads / FASTQ | `record(meta: Record, r1: Path?, r2: Path?, se: Path?, lr: Path?)` |
   | Assembly + reads | Assembly record + reads on separate lines |
   | Download / no input | No record input |
   | Alignment | `record(meta: Record, aln: Path)` |

   **b. Additional inputs** -- "Does this tool require a database or other additional files?"
   - Database: add `db: Path` as a separate input line
   - Optional files: add `optional_file: Path?`

   **c. Output files** -- "What output files does this tool produce?"
   - Ask for file extensions, descriptions, and whether they're single files or multiple
   - Use `file()` for single known files, `files()` for wildcards/multiple

   **d. Parameters** -- "Which command-line options should be exposed as Nextflow parameters?"
   - **Every parameter MUST be prefixed with the tool name**: `{tool}_{param}`
   - Ask for: name, type (string/integer/number/boolean), default value, description, CLI flag

   **e. Resource label:**
   | Label | When to use |
   |-------|-------------|
   | `process_low` | Default for most tools (4 CPU, 8GB, 4h) |
   | `process_medium` | BLAST-based, database searches (8 CPU, 32GB, 12h) |
   | `process_high` | Memory-intensive (12 CPU, 64GB, 24h) |
   | `process_single` | Single-threaded only (1 CPU, 4GB, 2h) |

   **f. Version command** -- "How does the tool report its version?"

   **g. Compressed input** -- "Does the tool handle gzipped (.gz) input natively?"

   **h. Citation key** and **keywords** for GroovyDoc

3. **Summarize everything** and ask the user to confirm before generating files.

---

### Phase 3: File Generation

**Goal:** Generate the 6 module files using `bactopia-scaffold`.

1. Construct the JSON config from the design decisions. Write it to `/tmp/scaffold-config.json`:

   ```json
   {
       "tool": "{tool_name}",
       "display_name": "{DisplayName}",
       "description": "{One-sentence description}",
       "process_name": "{TOOL_NAME}",
       "package": "{package_name}",
       "version": "{version}",
       "build": "{build}",
       "home_url": "{github_url}",
       "input_type": "assembly",
       "has_database": false,
       "handles_gz": false,
       "layout": "flat",
       "resource_label": "process_low",
       "version_command": "{version_command}",
       "citation_key": "{citation_key}",
       "keywords": ["{keyword1}", "{keyword2}"],
       "aggregation": {"strategy": "none"},
       "outputs": [
           {"name": "{field}", "extension": "{ext}", "description": "{desc}"}
       ],
       "parameters": [
           {"name": "{tool}_{param}", "type": "{type}", "default": "{default}", "description": "{desc}", "flag": "{--flag}"}
       ],
       "container_refs": {
           "toolName": "{from lookup}",
           "docker": "{from lookup}",
           "image": "{from lookup}"
       }
   }
   ```

2. Run the scaffold command:
   ```bash
   bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-scaffold.sh module --config /tmp/scaffold-config.json --bactopia-path . --pretty
   ```

3. The command creates 6 files:
   - `modules/{tool}/main.nf`
   - `modules/{tool}/module.config`
   - `modules/{tool}/schema.json`
   - `modules/{tool}/tests/main.nf.test`
   - `modules/{tool}/tests/nextflow.config`
   - `modules/{tool}/tests/nf-test.config`

---

### Phase 4: Review & Customize

**Goal:** Review generated files and make tool-specific adjustments.

1. **Module `main.nf`** -- the shell script block is a placeholder. Customize:
   - The actual tool command, flags, and I/O handling
   - Input decompression logic (if the tool doesn't handle .gz)
   - Database extraction logic (if database-dependent)
   - Any cleanup steps
   - Version extraction command

2. **Module `module.config`** -- review the `ext.args` construction:
   - Verify boolean/string/integer flag handling is correct for each parameter
   - Add any fixed flags (e.g., `--threads ${task.cpus}`)

3. **`schema.json`** -- verify parameter types and defaults match module.config

4. **Test files** -- verify test data paths and snapshot fields

5. **List all created files** with full paths.

6. **Remind the user** of follow-up steps:
   - If this module is part of a subworkflow, use `/add-subworkflow` next
   - If this is a standalone bactopia-tool, use `/add-bactopia-tool` instead (it handles all tiers)

---

## Edge Cases

1. **Multi-process modules**: For nested layouts (`modules/{tool}/run/`, `modules/{tool}/summary/`), the scaffold command currently generates flat layout. Manually move files to the nested structure after generation.

2. **No build string**: Container URLs will contain `TODO_BUILD` placeholders.

3. **Multi-package tools** (mulled containers): Container URLs cannot be auto-constructed. Flag for manual review.
