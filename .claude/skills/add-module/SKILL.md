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

3. Use WebFetch to query the Anaconda API:
   ```
   https://api.anaconda.org/package/bioconda/{package_name}
   ```
   If the package is not found on bioconda, try conda-forge:
   ```
   https://api.anaconda.org/package/conda-forge/{package_name}
   ```

4. From the API response, extract:
   - `latest_version` -- the version string
   - `summary` -- one-line tool description
   - `home` -- GitHub/documentation URL
   - From the `files` array, find the build string for the latest version:
     - Iterate through `files` looking for entries where `version == latest_version`
     - Prefer entries with `attrs.subdir == "linux-64"`
     - Fall back to `attrs.subdir == "noarch"` if no linux-64 build exists
     - Extract the build string from the filename: `{package}-{version}-{build}.tar.bz2` -- the `{build}` part (e.g., `hdfd78af_0`)

5. Construct container references:
   - `ext.toolName`: `"bioconda::{package}={version}"`
   - `ext.docker`: `"biocontainers/{package}:{version}--{build}"`
   - `ext.image`: `"https://depot.galaxyproject.org/singularity/{package}:{version}--{build}"`

6. **Present findings to the user** and ask them to confirm before proceeding:
   - Package name, version, build
   - Tool summary
   - GitHub URL
   - Constructed container URLs

7. **Check if the module already exists**: verify `modules/{tool}/` or `modules/{tool}/{process}/` does not already exist. If it does, warn the user.

---

### Phase 2: Information Gathering

**Goal:** Determine tool behavior, inputs, outputs, and parameters through web lookup + user interview.

1. **Fetch the tool's documentation** using WebFetch on the `home` URL from Phase 1.
   - If the URL is a GitHub repo, fetch the main page (includes README)
   - Extract: command-line options, input file types, output files, version command
   - If WebFetch fails, skip to step 2 and ask the user directly

2. **Ask the user to confirm or modify** the following (present what you found from docs):

   **a. Input type** -- "What does this tool take as input?"
   | User says | Input signature |
   |-----------|----------------|
   | Assembly / contigs / FASTA | `(_meta: Map, assembly: Path): Record` |
   | Reads / FASTQ | `(_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record` |
   | Assembly + reads | Assembly record + reads on separate lines |
   | Download / no input | No record input |
   | Alignment | `(_meta: Map, aligned_fa: Path): Record` |

   **b. Additional inputs** -- "Does this tool require a database or other additional files?"
   - Database: add `db: Path` as a separate input line
   - Optional files: add `optional_file: Path?` (null when not provided, use `!= null` checks in script block)

   **c. Output files** -- "What output files does this tool produce?"
   - Ask for file extensions, descriptions, and whether they're single files or multiple
   - Use `file()` for single known files, `files()` for wildcards/multiple

   **d. Parameters** -- "Which command-line options should be exposed as Nextflow parameters?"
   - **Every parameter MUST be prefixed with the tool name**: `{tool}_{param}` (e.g., `nohuman_db`, `mlst_scheme`, `bakta_skip_trna`)
   - This applies to ALL parameters: flags, paths, thresholds, database paths, etc.
   - The prefix prevents name collisions when multiple modules are loaded together
   - Ask for: name, type (string/integer/number/boolean), default value, description

   **e. Resource label** -- "How resource-intensive is this tool?"
   | Label | When to use |
   |-------|-------------|
   | `process_low` | Default for most tools (4 CPU, 8GB, 4h) |
   | `process_medium` | BLAST-based, database searches (8 CPU, 32GB, 12h) |
   | `process_high` | Memory-intensive (12 CPU, 64GB, 24h) |
   | `process_single` | Single-threaded only (1 CPU, 4GB, 2h) |

   **f. Version command** -- "How does this tool report its version?"
   - Common pattern: `tool --version 2>&1 | sed 's/pattern//'`
   - If no CLI version support: use hardcoded VERSION + `ext.version` in module.config

   **g. Compressed input** -- "Does this tool accept gzipped (.gz) input files natively?"
   - If no: include gzip decompression block in the script
   - If yes: skip the decompression block

   **h. Process name** -- Derive the UPPER_CASE process name:
   - Simple module: `TOOL_NAME` (e.g., `MLST`, `SISTR`)
   - Multi-process: `TOOL_PROCESS` (e.g., `BAKTA_RUN`, `BAKTA_DOWNLOAD`)

   **i. Citation key** -- Ask for the citation key for `@citation` tag

3. **Summarize all gathered information** and ask the user to confirm before generating files.

---

### Phase 3: File Generation

**Goal:** Generate all module files using the templates below.

Create the directory structure first:
- Simple module: `modules/{tool}/` and `modules/{tool}/tests/`
- Multi-process: `modules/{tool}/{process}/` and `modules/{tool}/{process}/tests/`

Generate these files:

#### File 1: main.nf

Use the appropriate template based on input type.

**Template A: Assembly input (most common)**

```groovy
/**
 * {One-sentence summary}.
 *
 * Uses [{ToolName}]({github_url}) to {detailed description of what the tool does
 * and the outputs it produces}.
 *
 * @status stable
 * @keywords {comma, separated, keywords}
 * @tags complexity:{simple|moderate|complex} input-type:{none|single|multiple} output-type:{single|multiple} features:{feature1,feature2}
 * @citation {citation_key}
 *
 * @note {Optional: Database Required, Database Bundled, etc.}
 * {Optional: Additional context about requirements.}
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, {field1}, {field2}, results, logs, nf_logs, versions)
 * - `{field1}`: {Description of tool-specific output}
 * - `{field2}`: {Description of tool-specific output}
 */
nextflow.preview.types = true

process {PROCESS_NAME} {
    tag "${prefix}"
    label '{process_label}'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        {field1}: file("${prefix}.{ext1}"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.{ext1}")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    {tool_command} \\
        ${task.ext.args} \\
        {input_flag} ${assembly_name} \\
        {output_flags}

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        {tool}: $( {version_command} )
    END_VERSIONS
    """
}
```

**Template B: Assembly + database input**

Same as Template A but add database input:
```groovy
    input:
    (_meta: Map, assembly: Path): Record
    db: Path
```

And add database handling in the script block:
```groovy
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        DB_PATH=database
    else
        DB_PATH=${db}
    fi

    {tool_command} \\
        ${task.ext.args} \\
        --db $DB_PATH \\
        ...
    """
```

**Template C: Reads input**

```groovy
    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record

    // ... in script block:
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    meta.single_end = has_se && !has_r1 && !has_r2

    read_inputs = meta.single_end ? "${se}" : "${r1} ${r2}"
```

**Template D: Download module (no sample input)**

```groovy
    // No record input -- download modules fetch from external sources

    output:
    record(
        db: files("{output_dir}/*"),
        logs: files("*.{log,err}", optional: true)
    )
```

**Important rules for main.nf:**
- ALWAYS include `nextflow.preview.types = true` before the process
- ALWAYS construct a fresh `meta` map in the script block (never modify `_meta`)
- ALWAYS emit `record()` with named fields + results/logs/nf_logs/versions
- Use `file()` for single known files, `files()` for wildcards/multiple files
- Use 4 spaces for indentation
- If the tool has no `--version` flag, use hardcoded VERSION:
  ```groovy
  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
  def VERSION = '{version}'
  ```
- Skip the gzip decompression block if the tool handles .gz natively
- Use exactly **one space** after the colon in record fields (e.g., `meta: meta`, not `meta:  meta`)
- Always include `# Cleanup` comment before the versions block (even if no cleanup needed)

#### File 2: module.config

```groovy
params {
    // {tool}
    {tool}_{param1} = {default1}
    {tool}_{param2} = {default2}
}

process {
    withName: '{PROCESS_NAME}' {
        ext.wf = params.wf
        ext.scope = "sample"
        ext.subdir = ""
        ext.logs_subdir = ""
        ext.process_name = "{tool}"

        // Tool arguments
        ext.args = [
            params.{tool}_{param1} ? "--flag1 ${params.{tool}_{param1}}" : "",
            "--param2 ${params.{tool}_{param2}}"
        ].join(' ').replaceAll("\\s{2,}", " ").trim()

        // Environment information
        ext.toolName = "bioconda::{package}={version}".replace("=", "-").replace(":", "-").replace(" ", "-")
        ext.docker = "biocontainers/{package}:{version}--{build}"
        ext.image = "https://depot.galaxyproject.org/singularity/{package}:{version}--{build}"
        ext.condaDir = "${params.condadir}"
    }
}
```

**Block ordering:** routing first, then `// Tool arguments`, then `// Environment information`, then `// Module-specific parameters` (optional).

**Params conventions:** parameters in alphabetical order, section comment uses snake_case (`// {tool}`), capitalize `// No parameters` for empty blocks.

If no parameters:
```groovy
params {
    // No parameters
}

process {
    withName: '{PROCESS_NAME}' {
        ext.wf = params.wf
        ext.scope = "sample"
        ext.subdir = ""
        ext.logs_subdir = ""
        ext.process_name = "{tool}"

        // Tool arguments
        ext.args = ""

        // Environment information
        // ... same as above
    }
}
```

If no CLI version command, add:
```groovy
        // Version information not provided by tool on CLI
        ext.version = "{version}"
```

#### File 3: schema.json

```json
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://raw.githubusercontent.com/bactopia/bactopia/master/modules/{module_path}/schema.json",
    "title": "{Tool Name} Module",
    "description": "A module for {brief description}",
    "type": "object",
    "$defs": {
        "{tool}_parameters": {
            "title": "{Tool Name} Parameters",
            "type": "object",
            "description": "",
            "default": "",
            "fa_icon": "fas fa-exclamation-circle",
            "properties": {
                "{tool}_{param1}": {
                    "type": "{string|integer|number|boolean}",
                    "default": {default_value},
                    "description": "{Parameter description}",
                    "fa_icon": "{icon_by_type}"
                }
            }
        }
    },
    "allOf": [
        {
            "$ref": "#/$defs/{tool}_parameters"
        }
    ]
}
```

- `{module_path}` is the relative path under `modules/` (e.g., `mlst` or `bakta/run`)
- If no parameters, the `properties` object should be empty `{}`
- `{icon_by_type}`: `string` = `fas fa-font`, `integer` = `fas fa-hashtag`, `number` = `fas fa-percentage`, `boolean` = `fas fa-toggle-on`

#### File 4: tests/main.nf.test

**Assembly input test:**

```groovy
nextflow_process {
    name "Test {PROCESS_NAME}"
    script "../main.nf"
    process "{PROCESS_NAME}"
    tag "modules"
    tag "{tool}"

    test("{tool} - module - GCF_000292685") {
        when {
            params {
                test_data_dir = System.getenv("BACTOPIA_TESTS") ?: ""
            }
            process {
                """
                input[0] = Channel.of(
                    record(
                        _meta: [name: "GCF_000292685"],
                        assembly: file("${params.test_data_dir}/data/species/portiera/genome/GCF_000292685.fna.gz")
                    )
                )
                """
            }
        }

        then {
            def record = process.out[0][0]
            assertAll(
                { assert process.success },
                { assert snapshot(
                    record.meta,
                    record.{output_field},
                    record.versions
                ).match() }
            )
        }
    }
}
```

If the tool has a database input, add it:
```groovy
                input[1] = file("${params.test_data_dir}/data/datasets/{tool}/{db_file}")
```

**Reads input test:**
```groovy
                input[0] = Channel.of(
                    record(
                        _meta: [name: "SRR2838702"],
                        r1: file("${params.test_data_dir}/data/species/portiera/illumina/SRR2838702_R1.fastq.gz"),
                        r2: file("${params.test_data_dir}/data/species/portiera/illumina/SRR2838702_R2.fastq.gz"),
                        se: null,
                        lr: null
                    )
                )
```

**Snapshot rules:**
- ALWAYS include: `record.meta`, tool-specific output fields, `record.versions`
- NEVER include: `record.results`, `record.logs`, `record.nf_logs` (unstable across runs)

#### File 5: tests/nextflow.config

```groovy
// Minimal config for module-level testing of {PROCESS_NAME}
nextflow.preview.types = true
nextflow.enable.strict = true

params {
    workflow {
        name = "{tool}"
        logo_name = "bactopia-tools"
        description = "{Tool description}"
        ext = "fna"
    }

    bactopia_version = '4.0.0'
    bactopia_cache = System.getenv("BACTOPIA_CACHEDIR") ?: "${System.getenv('HOME')}/.bactopia"
    condadir = "${params.bactopia_cache}/conda"
    wf = params.workflow.name
    merge_folder = "merged-results"
    test_data_dir = System.getenv("BACTOPIA_TESTS") ?: ""
    is_ci = true

    // Max Job Request Parameters
    max_retry = 1
    max_time = 2.h
    max_memory = 8.GB
    max_cpus = 2

    // Nextflow Profile Parameters
    registry = "quay.io"
    singularity_cache = "${params.bactopia_cache}/singularity"
    singularity_pull_docker_container = false
    container_opts = ""
}

includeConfig "../module.config"
includeConfig "{depth}conf/base.config"
includeConfig "{depth}conf/profiles.config"
```

**Path depth:**
- Simple module (`modules/{tool}/tests/`): `{depth}` = `../../../`
- Multi-process (`modules/{tool}/{process}/tests/`): `{depth}` = `../../../../`

#### File 6: tests/nf-test.config

This is identical across all modules:

```groovy
config {
    testsDir "."
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"
    configFile "nextflow.config"
    profile "docker"
    options "--is_ci --max_memory 8.GB"

    plugins {
        load "nft-utils@0.0.5"
    }
}
```

---

### Phase 4: Summary and Next Steps

After generating all files, present:

1. **List of created files** with full paths
2. **Remaining integration steps** (NOT done by this skill):
   - Add `includeConfig "./modules/{module_path}/module.config"` to `nextflow.config`
   - Add citation to `data/citations.yml` if needed
   - Create subworkflow (use future `add-subworkflow` skill)
   - Create workflow entry point (use future `add-workflow` skill)
   - Run tests: `cd modules/{module_path}/tests && nf-test test main.nf.test --update-snapshot`

---

## Edge Cases

1. **Package not found on bioconda**: Try conda-forge. If still not found, ask user for version/build manually and note that containers may not be available on biocontainers.

2. **No build string in API response**: Warn user and create module.config with a `TODO` placeholder for the build hash. The module will not work until correct container URLs are provided.

3. **No --version CLI support**: Use the hardcoded VERSION pattern (see Template A notes). Add `ext.version` to module.config.

4. **Multi-package tools** (e.g., mulled containers): Warn the user that container URLs cannot be auto-constructed from a single package name. Flag for manual review -- the user will need to find the correct mulled container hash.

5. **Module directory already exists**: Warn the user before creating files. Ask if they want to overwrite or create in a subdirectory.

6. **WebFetch failures**: If Anaconda API or GitHub docs cannot be fetched, fall back to asking the user for all information manually. Do not block on WebFetch failures.

## Test Data Reference

| Species | Path | When to use |
|---------|------|-------------|
| **Portiera** | `portiera/genome/GCF_000292685.fna.gz` | **Default for all assembly modules** |
| S. aureus | `staphylococcus_aureus/genome/GCF_000017085.fna` | Species-specific typing tools |
| N. gonorrhoeae | `neisseria_gonorrhoeae/genome/GCF_001047255.fna` | Has .gz variant |
| H. influenzae | `haemophilus_influenzae/genome/GCF_900478275.fna` | Has .gz variant |
| Portiera reads | `portiera/illumina/SRR2838702_R{1,2}.fastq.gz` | **Default for read-based modules** |

All paths are relative to `$BACTOPIA_TESTS/data/species/`.
