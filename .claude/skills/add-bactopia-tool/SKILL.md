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

## Phased Workflow

Follow these phases in order. When unsure about ANYTHING, ask the user rather than guess.

---

### Phase 1: Package Verification

**Goal:** Confirm the bioconda package exists and retrieve version/container information.

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

---

### Phase 2: Tool Design

**Goal:** Gather all design decisions in one pass so files can be generated coherently.

1. **Fetch the tool's documentation** using WebFetch on the `home` URL from Phase 1.
   - Extract: command-line options, input file types, output files, version command
   - If WebFetch fails, ask the user directly

2. **Ask the user to confirm or modify** the following (present what you found from docs):

   **a. Tool identity:**
   - Tool name (snake_case, used for directories and param prefixes)
   - Display name (for GroovyDoc headers)
   - One-sentence description

   **b. Input type** -- determines which BACTOPIATOOL_INIT channel to use:

   | Input Type | Channel | `params.workflow.ext` | Module record input |
   |---|---|---|---|
   | Assembly | `assembly` | `['fna']` | `record(meta: Record, fna: Path)` |
   | Reads | `reads` | `['fastq']` | `record(meta: Record, r1: Path?, r2: Path?, se: Path?, lr: Path?)` |
   | Assembly + reads | `assembly_reads` | `['fna', 'fastq']` | `record(meta: Record, fna: Path, r1: Path?, r2: Path?, se: Path?, lr: Path?)` |
   | Proteins | `proteins` | `['faa']` | `record(meta: Record, faa: Path)` |
   | GFF | `gff` | `['gff']` | `record(meta: Record, gff: Path)` |
   | GenBank | `gbff` | `['gbk']` | `record(meta: Record, gbff: Path)` |

   **c. Additional inputs** -- "Does this tool require a database or other files?"
   - Database: add `db: Path` as a separate input (affects subworkflow take + workflow params)
   - Other files: add as typed inputs

   **d. Output files** -- "What output files does this tool produce?"
   - File extensions, descriptions, single vs multiple
   - Which field to use for aggregation (e.g., `tsv`, `csv`, `report`)

   **e. Parameters** -- "Which command-line options should be exposed?"
   - Every parameter MUST be prefixed with the tool name: `{tool}_{param}`
   - Ask for: name, type (string/integer/number/boolean), default value, description, CLI flag

   **f. Resource label:**
   | Label | When to use |
   |-------|-------------|
   | `process_low` | Default for most tools (4 CPU, 8GB, 4h) |
   | `process_medium` | BLAST-based, database searches (8 CPU, 32GB, 12h) |
   | `process_high` | Memory-intensive (12 CPU, 64GB, 24h) |
   | `process_single` | Single-threaded only (1 CPU, 4GB, 2h) |

   **g. Version command** -- how does this tool report its version?

   **h. Compressed input** -- does this tool accept gzipped (.gz) input natively?

   **i. Aggregation strategy:**
   - **CSVTK_CONCAT** (default, most common) -- which output field, tsv or csv format
   - **Dedicated summary module** (rare -- only when tool has its own aggregation command)
   - **No aggregation** -- tool doesn't produce per-sample tabular output

   **j. Citation info** -- citation key, tool name, URL, description, citation text

   **k. Keywords** for GroovyDoc

   **l. Test data** -- which species/sample to use:
   - Default: `portiera/GCF_000292685` for assembly tools
   - Species-specific typing tools: `staphylococcus_aureus/GCF_000017085`
   - Database-dependent: which test dataset path

3. **Summarize everything** and ask the user to confirm before generating files.

---

### Phase 3: File Generation

**Goal:** Generate all 15 files across the three tiers using `bactopia-scaffold`.

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

3. The command creates all 15 files. Review the output to confirm which files were created.

---

### Phase 4: Review & Customize

**Goal:** Review generated files and make tool-specific adjustments.

The templates produce correct scaffolds but many tools need customization:

1. **Module `main.nf`** -- the shell script block is a placeholder. Customize:
   - The actual tool command, flags, and I/O handling
   - Input decompression logic (if the tool doesn't handle .gz)
   - Database extraction logic (if database-dependent)
   - Any cleanup steps
   - Version extraction command

2. **Module `module.config`** -- review the `ext.args` construction:
   - Verify boolean/string/integer flag handling is correct for each parameter
   - Add any fixed flags (e.g., `--threads ${task.cpus}`)

3. **Subworkflow `main.nf`** -- usually correct as-is for CSVTK_CONCAT pattern. Check:
   - The `@input` GroovyDoc matches the subworkflow's input name (may differ from module input)
   - The `@output` field descriptions are accurate

4. **Workflow `main.nf`** -- check GroovyDoc `@publish` sections match actual outputs

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

3. **Remind the user** to run these follow-up skills:
   - `/update-catalog` -- regenerate `catalog.json` and `llms.txt`
   - `/merge-schemas` on the new workflow -- generate `nextflow_schema.json`
   - `/run-tests` -- verify the scaffolded tests pass

4. **Note:** `nextflow_schema.json` is NOT generated by this skill -- `/merge-schemas` handles it automatically from the module `schema.json` files.

---

## Edge Cases

1. **Package not found**: The lookup command tries bioconda first, then conda-forge. If both fail, ask the user for version/build manually.

2. **No build string**: Container URLs will contain `TODO_BUILD` placeholders. Flag for manual review.

3. **No --version CLI support**: Use hardcoded VERSION pattern in the module main.nf.

4. **Multi-package tools** (mulled containers): Warn the user that container URLs cannot be auto-constructed. Flag for manual review.

5. **Component already exists**: The lookup output includes `existing_components`. Warn before proceeding.

## Test Data Reference

| Species | Path | When to use |
|---------|------|-------------|
| **Portiera** | `portiera/compressed/GCF_000292685/main/assembler/GCF_000292685.fna.gz` | **Default for assembly tools** |
| S. aureus | `staphylococcus_aureus/compressed/GCF_000017085/main/assembler/GCF_000017085.fna.gz` | Species-specific typing tools |
| Portiera reads | `portiera/illumina/SRR2838702_R{1,2}.fastq.gz` | **Default for read-based tools** |

All paths are relative to `$BACTOPIA_TESTS/species/`.
