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

## Interactive Questioning

This skill is interactive -- ask the user early and often, especially before creating files.

- **Multiple questions at once:** Use `AskUserQuestion` popups (up to 4 questions per batch).
  Mark the recommended option with "(Recommended)" at the end of its label and place it first.
- **Single simple question:** Just ask in chat, no popup needed.
- **When in doubt:** Ask. It's cheaper to clarify upfront than to regenerate files.

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

### Phase 2: Tool Design

**Goal:** Gather all design decisions using interactive prompts so files can be generated coherently.

**Important:** Use the `AskUserQuestion` tool for structured choices throughout this phase.
Present up to 4 questions per batch. Mark the recommended option (based on WebFetch findings)
with "(Recommended)" at the end of its label and place it first in the options list.

1. **Fetch the tool's documentation** using WebFetch on the `home` URL from Phase 1.
   - Extract: command-line options, input file types, output files, version command
   - If WebFetch fails, ask the user directly

2. **Batch 1: Core design choices** (AskUserQuestion, up to 4 questions)

   Based on WebFetch findings, ask these structured questions:

   **Question 1 -- Input type:**
   | Input Type | Record fields |
   |---|---|
   | Assembly | `record(meta: Record, fna: Path)` |
   | Reads | `record(meta: Record, r1: Path?, r2: Path?, se: Path?, lr: Path?)` |
   | Assembly + reads | Assembly record + reads on separate lines |
   | Alignment | `record(meta: Record, aln: Path)` |
   | Download / no input | No record input |

   Options (pick top 3 most relevant, "Other" is auto-added for the rest):
   - Assembly -- takes FASTA assembly files
   - Reads -- takes FASTQ read files
   - Assembly + Reads -- takes both FASTA and FASTQ

   **Question 2 -- Database requirement:**
   - No database needed
   - Yes, requires a user-provided database

   **Question 3 -- Resource label:**
   - process_low -- 4 CPU, 8GB, 4h (default for most tools)
   - process_medium -- 8 CPU, 32GB, 12h (BLAST-based, database searches)
   - process_high -- 12 CPU, 64GB, 24h (memory-intensive)
   - process_single -- 1 CPU, 4GB, 2h (single-threaded only)

   **Question 4 -- Compressed input:**
   - Yes, handles .gz natively
   - No, needs decompression first

3. **Run test-data discovery** based on the input type selected in Batch 1:
   ```bash
   bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-scaffold.sh test-data --input-type {input_type} --bactopia-path . --pretty
   ```
   This returns species/accession combinations already used by similar modules, with
   pre-computed `test_data_path`, `test_uncompressed_path`, `test_species`, and
   `test_sample_id` values. Use the returned paths directly in the scaffold config
   (Phase 3) -- do NOT construct paths manually.

4. **Batch 2: Outputs, parameters, and test data** (AskUserQuestion, up to 3 questions)

   **Question 1 -- Output files:**
   Present what WebFetch found. Ask the user to confirm file extensions, descriptions,
   and whether each is single (`file()`) or multiple (`files()`).

   **Question 2 -- User parameters:**
   Only flags representing user-meaningful analysis choices (identity thresholds, scheme
   selection, algorithm toggles). Exclude infrastructure params (see below).

   **Question 3 -- Test data species:**
   Present top 3-4 species from the test-data discovery results. Recommend species that
   exercise the tool's functionality. Include the accession in each option's description.

5. **Present auto-detected details for confirmation.**

   After the structured choices, present these findings from WebFetch in a summary and
   ask the user to confirm or request changes:

   - **Tool identity**: name (snake_case), display name, one-sentence description
   - **Version command**: how the tool reports its version
   - **Citation key** and **keywords** for GroovyDoc

   **Infrastructure vs. user parameters (do NOT expose these):**

   | Tool flag | Wired to | Where |
   |-----------|----------|-------|
   | `--prefix`, `--label`, `--sample-name`, etc. | `prefix` variable (`task.ext.prefix ?: "${_meta.name}"`) | Shell block |
   | `--threads`, `--cpus`, `-t`, `-p`, etc. | `${task.cpus}` | Shell block or `ext.args` in module.config |
   | `--output`, `--outdir`, `-o`, etc. | Usually `.` or `${prefix}` | Shell block |

   These are written directly in the module's shell block (e.g., `--prefix ${prefix}`,
   `--threads ${task.cpus}`). The `prefix` variable is set in every module's script block
   as `prefix = task.ext.prefix ?: "${_meta.name}"` and carries the sample name.

   Only expose flags that represent **user-meaningful analysis choices**.

   Every user parameter MUST be prefixed with the tool name: `{tool}_{param}`.

   **Parameter defaults:**
   - Do NOT assume a default is needed. Ask the user whether each parameter should have
     a specific default value.
   - If a string parameter needs a default, use an empty string `""`, never `null`.
   - Only include a parameter in the config's `"parameters"` array if the user confirms
     it should be exposed.

6. **Final confirmation** (AskUserQuestion, 1 question)

   After presenting the summary, ask:
   - Looks good, proceed to file generation
   - I need to make changes (user provides details via "Other" or notes)

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
       },
       "test_species": "{species}",
       "test_sample_id": "{sample_id}",
       "test_data_path": "{compressed_path}",
       "test_uncompressed_path": "{uncompressed_path}"
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
   - **Always preserve the `# Cleanup` comment line** -- even if empty, it marks where
     cleanup steps go and keeps the shell block structure consistent across all modules
   - Version extraction command

2. **Module `module.config`** -- review the `ext.args` construction:
   - Verify boolean/string/integer flag handling is correct for each parameter
   - Add any fixed flags (e.g., `--threads ${task.cpus}`)

3. **`schema.json`** -- verify parameter types and defaults match module.config

4. **Test files** -- verify test data paths and snapshot fields

5. **Run the linter** to catch structural issues before proceeding:
   ```bash
   bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-lint.sh {tool} --bactopia-path .
   ```
   This runs `bactopia-lint` scoped to the new module. Fix any FAILs before moving on.
   Common issues:
   - JS005: type/default mismatch in schema.json (e.g., `type=string` but `default=null`)
   - M035: citation key not found in `data/citations.yml`

6. **Update `data/citations.yml`** -- add the tool citation entry in alphabetical order:
   ```yaml
   {tool}:
     name: "{ToolName}"
     link: "{github_url}"
     description: "{One-sentence description}"
     cite: "{Full citation text}"
   ```

7. **List all created files** with full paths.

8. **Remind the user** to run these follow-up steps in order:
   1. `/run-tests {tool} module --generate` -- generate snapshots and verify the module test passes (new modules have no existing snapshots)
   2. If this module is part of a subworkflow, use `/add-subworkflow` next
   3. If this is a standalone bactopia-tool, use `/add-bactopia-tool` instead (it handles all tiers)

   The `--generate` flag is required because newly scaffolded modules have no
   snapshot files yet. Without it, nf-test will fail immediately on missing
   snapshots.

---

## Edge Cases

1. **Multi-process modules**: For nested layouts (`modules/{tool}/run/`, `modules/{tool}/summary/`), the scaffold command currently generates flat layout. Manually move files to the nested structure after generation.

2. **No build string**: Container URLs will contain `TODO_BUILD` placeholders.

3. **Multi-package tools** (mulled containers): Container URLs cannot be auto-constructed. Flag for manual review.

## Test Data Discovery

Test data paths are discovered dynamically from existing module tests using:
```bash
bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-scaffold.sh test-data --input-type {type} --bactopia-path . --pretty
```

This scans `modules/*/tests/main.nf.test` for paths matching the input type and returns
pre-computed template variables. Always use the discovered paths -- never construct test
data paths manually. The output includes `test_data_path` (compressed, for subworkflow
tests), `test_uncompressed_path` (for module tests), `test_species`, and `test_sample_id`.

Supported input types: `assembly`, `reads`, `assembly_reads`, `proteins`, `gff`, `genbank`.
