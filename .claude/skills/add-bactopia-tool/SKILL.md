---
name: add-bactopia-tool
description: >-
  Scaffold a complete Bactopia Tool across all three tiers -- module, subworkflow, and
  workflow entry point under workflows/bactopia-tools/. Creates all files (main.nf,
  module.config, schema.json, nextflow.config, tests) for the common single-tool pattern.
  Use when asked to add a new bactopia tool, create a bactopia tool, scaffold a complete
  tool, add a new analysis tool to bactopia-tools, or wire up a bioconda package as a
  bactopia-tool. This skill handles the full pipeline from package lookup through file
  generation -- do not use add-module or add-subworkflow separately when the goal is a
  complete bactopia-tool.
---

# Add Bactopia Tool

Scaffold a complete Bactopia Tool pipeline from a bioconda/conda-forge package. This creates **all three tiers** in one shot:

1. **Module** (`modules/{tool}/`) -- the Nextflow process that runs the tool
2. **Subworkflow** (`subworkflows/{tool}/`) -- orchestrates the module + aggregation
3. **Workflow** (`workflows/bactopia-tools/{tool}/`) -- user-facing entry point

This skill handles the **common single-tool pattern** which covers ~80% of bactopia-tools (abricate, mlst, bakta, quast, sistr, etc.). Multi-stage pipelines like snippy and pangenome should be hand-built.

## Prerequisites

Before using this skill, read:
- `.claude/docs/standards/05-module-documentation.md` -- Module GroovyDoc standards
- `.claude/docs/standards/04-subworkflow-documentation.md` -- Subworkflow GroovyDoc standards

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

**Goal:** Confirm the package exists on bioconda and retrieve version/container information.

1. Ask the user for the **bioconda package name** (e.g., `mlst`, `bakta`, `ssuissero`).

2. Run the lookup command:
   ```bash
   bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-scaffold.sh lookup {package_name} --bactopia-path . --pretty
   ```

3. The output includes:
   - `package`, `channel`, `version`, `build` -- package identity
   - `summary`, `home` -- tool description and documentation URL
   - `container_refs` -- `toolName`, `docker`, `image` strings
   - `existing_components` -- which of module/subworkflow/workflow already exist

4. **Present findings to the user** and ask them to confirm before proceeding.
   If any existing components are found, warn the user.
   If the package is not found, ask the user to verify the name.
   If the package does not exist on bioconda, inform the user you cannot proceed until a valid bioconda package is provided.

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
   Determines which BACTOPIATOOL_INIT channel to use.

   | Input Type | Channel | `params.workflow.ext` | Module record input |
   |---|---|---|---|
   | Assembly | `assembly` | `['fna']` | `record(meta: Record, fna: Path)` |
   | Reads | `reads` | `['fastq']` | `record(meta: Record, r1: Path?, r2: Path?, se: Path?, lr: Path?)` |
   | Assembly + reads | `assembly_reads` | `['fna', 'fastq']` | `record(meta: Record, fna: Path, r1: Path?, r2: Path?, se: Path?, lr: Path?)` |
   | Proteins | `proteins` | `['faa']` | `record(meta: Record, faa: Path)` |
   | GFF | `gff` | `['gff']` | `record(meta: Record, gff: Path)` |
   | GenBank | `gbff` | `['gbk']` | `record(meta: Record, gbff: Path)` |

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

4. **Batch 2: Aggregation and test data** (AskUserQuestion, up to 2 questions)

   **Question 1 -- Aggregation strategy:**
   - CSVTK_CONCAT -- concatenate per-sample tabular output (most common)
   - Dedicated summary module -- tool has its own aggregation command (rare)
   - No aggregation -- tool doesn't produce per-sample tabular output

   **Question 2 -- Test data species:**
   Present top 3-4 species from the test-data discovery results. Recommend species that
   exercise the tool's functionality (e.g., species with multiple MLST schemes for
   typing tools, species with known resistance genes for AMR tools). Include the accession
   in each option's description.

5. **Present auto-detected details for confirmation.**

   After the structured choices, present these findings from WebFetch in a summary and
   ask the user to confirm or request changes:

   - **Tool identity**: name (snake_case), display name, one-sentence description
   - **Output files**: extensions, descriptions, aggregation field
   - **User parameters**: only flags representing user-meaningful analysis choices
     (identity thresholds, scheme selection, algorithm toggles). Exclude infrastructure
     params (see below).
   - **Version command**: how the tool reports its version
   - **Citation**: key, name, URL, description, citation text
   - **Keywords**: for GroovyDoc

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

**Goal:** Generate all 16 files across the three tiers using `bactopia-scaffold`.

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
       "input_type": "{assembly|reads|assembly_reads|proteins|gff|genbank}",
       "has_database": false,
       "handles_gz": false,
       "layout": "flat",
       "resource_label": "{process_low|process_medium|process_high|process_single}",
       "version_command": "{version_command}",
       "citation_key": "{citation_key}",
       "keywords": ["{keyword1}", "{keyword2}"],
       "aggregation": {
           "strategy": "{csvtk_concat|dedicated_summary|none}",
           "field": "{output_field}",
           "format": "{tsv|csv}"
       },
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
       "test_uncompressed_path": "{uncompressed_path}",
       "test_dataset": "{dataset_path_or_empty}",
       "test_dataset2": "",
       "test_dataset3": ""
   }
   ```

   For database-dependent tools, also include:
   ```json
   {
       "database": {
           "param_name": "{tool}_db",
           "test_path": "datasets/{tool}/{db_file}"
       }
   }
   ```

2. Run the scaffold command:
   ```bash
   bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-scaffold.sh tool --config /tmp/scaffold-config.json --bactopia-path . --pretty
   ```

3. The command creates all 16 files. Review the output to confirm which files were created.

---

### Phase 4: Review & Customize

**Goal:** Review generated files and make tool-specific adjustments.

The templates produce correct scaffolds but many tools need customization:

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

3. **Subworkflow `main.nf`** -- usually correct as-is for CSVTK_CONCAT pattern. Check:
   - The `@input` GroovyDoc matches the subworkflow's input name (may differ from module input)
   - The `@output` field descriptions are accurate

4. **Workflow `main.nf`** -- check GroovyDoc `@publish` sections match actual outputs

5. **Run the linter** to catch structural issues before proceeding:
   ```bash
   bash .claude/skills/add-bactopia-tool/scripts/run-bactopia-lint.sh {tool} --bactopia-path .
   ```
   This runs `bactopia-lint` scoped to the new module, subworkflow, and workflow.
   Fix any FAILs before moving to Phase 5. Common issues:
   - JS005: type/default mismatch in schema.json (e.g., `type=string` but `default=null`)
   - S011: misaligned include braces in subworkflow
   - M035/S019: citation key not found in `data/citations.yml`

---

### Phase 5: Integration & Next Steps

**Goal:** Wire up citations and inform the user about remaining steps.

1. **Update `data/citations.yml`** -- add the tool citation entry in alphabetical order:
   ```yaml
   {tool}:
     name: "{ToolName}"
     link: "{github_url}"
     description: "{One-sentence description}"
     cite: "{Full citation text}"
   ```

2. **List all created files** with full paths.

3. **Remind the user** to run these follow-up skills in order:
   1. `/run-tests {tool} module and subworkflow --generate` -- generate snapshots and verify tests pass (new tools have no existing snapshots)
   2. `/update-catalog` -- regenerate `catalog.json` and `llms.txt` (only after tests pass)
   3. `/merge-schemas` on the new workflow -- generate `nextflow_schema.json`
   4. `/run-tests {tool} workflow --generate` -- generate snapshots and verify the workflow test passes

   The `--generate` flag is required because newly scaffolded tools have no
   snapshot files yet. Without it, nf-test will fail immediately on missing
   snapshots.

   There is no point running `/update-catalog` or `/merge-schemas` if the
   module/subworkflow tests are failing.

4. **Note:** `nextflow_schema.json` is NOT generated by this skill -- `/merge-schemas` handles it automatically from the module `schema.json` files.

---

## Edge Cases

1. **Package not found**: The lookup command tries bioconda first, then conda-forge. If both fail, ask the user for version/build manually.

2. **No build string**: Container URLs will contain `TODO_BUILD` placeholders. Flag for manual review.

3. **No --version CLI support**: Use hardcoded VERSION pattern in the module main.nf.

4. **Multi-package tools** (mulled containers): Warn the user that container URLs cannot be auto-constructed. Flag for manual review.

5. **Component already exists**: The lookup output includes `existing_components`. Warn before proceeding.

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
