# task.ext Properties Reference

## Overview

Bactopia modules use `task.ext` properties to configure process behavior, container settings, and output organization. These properties are defined in `module.config` files and accessed within process scripts.

## Core Properties

### Container and Environment

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| `toolName` | String | Tool identifier for conda environment naming | `"bioconda-abricate-1.0.1"` |
| `docker` | String | Docker container reference | `"biocontainers/abricate:1.0.1--ha8f3691_1"` |
| `image` | String | Singularity/Apptainer image URL | `"https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1"` |
| `condaDir` | String | Conda environment directory | `"${params.condadir}"` |

### Output Organization

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| `scope` | String | Output scope: `"sample"` or `"run"` | `"sample"` |
| `process_name` | String | Process name for output paths | `"abricate"` |
| `subdir` | String | Subdirectory within output path | `"${params.abricate_db}"` |
| `logs_subdir` | String | Subdirectory for log files | `"abricate-run"` |
| `prefix` | String | Output filename prefix override | `"${meta.id}"` |

### Workflow Context

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| `wf` | String | Workflow identifier | `"bactopia"` or `"teton"` |

### Tool Arguments

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| `args` | String | Primary tool arguments | `"--db ncbi --minid 80"` |
| `args2` | String | Secondary tool arguments | `"-l S -t 0"` |

### Module-Specific Properties

Some modules define custom properties for internal logic:

| Property | Module | Description |
|----------|--------|-------------|
| `bracken_read_length` | bracken | Read length for Bracken estimation |
| `kraken2_keep_raw_output` | bracken | Keep raw Kraken2 output files |
| `kraken2_keep_filtered_reads` | bracken | Keep filtered read files |
| `skip_krona` | bracken | Skip Krona visualization |
| `max_retry` | various | Maximum retry attempts for failed processes |

## Configuration Patterns

### Basic Module Configuration

```groovy
process {
    withName: 'MLST' {
        // Tool arguments
        ext.args = [
            params.mlst_nopath ? "--nopath" : "",
            params.mlst_scheme ? "--scheme ${params.mlst_scheme}" : "",
            "--minid ${params.mlst_minid}",
            "--mincov ${params.mlst_mincov}"
        ].join(' ').replaceAll("\\s{2,}", " ").trim()

        // Container configuration
        ext.toolName = "bioconda::mlst=2.23.0".replace("=", "-").replace(":", "-").replace(" ", "-")
        ext.docker = "biocontainers/mlst:2.23.0--hdfd78af_1"
        ext.image = "https://depot.galaxyproject.org/singularity/mlst:2.23.0--hdfd78af_1"
        ext.condaDir = "${params.condadir}"

        // Output organization
        ext.wf = params.wf
        ext.scope = "sample"
        ext.subdir = ""
        ext.logs_subdir = ""
        ext.process_name = "mlst"
    }
}
```

### Multi-Process Configuration

When a tool has multiple processes (e.g., run + summary), use shared and process-specific blocks:

```groovy
process {
    // Shared configuration for all ABRICATE processes
    withName: 'ABRICATE_RUN|ABRICATE_SUMMARY' {
        ext.args = [
            "--db ${params.abricate_db}",
            "--minid ${params.abricate_minid}",
            "--mincov ${params.abricate_mincov}"
        ].join(' ').replaceAll("\\s{2,}", " ").trim()

        ext.toolName = "bioconda::abricate=1.0.1".replace("=", "-").replace(":", "-").replace(" ", "-")
        ext.docker = "biocontainers/abricate:1.0.1--ha8f3691_1"
        ext.image = "https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1"
        ext.condaDir = "${params.condadir}"

        ext.wf = params.wf
        ext.scope = "sample"
    }

    // Process-specific: ABRICATE_RUN
    withName: 'ABRICATE_RUN' {
        ext.subdir = params.abricate_db
        ext.logs_subdir = ""
        ext.process_name = "abricate"
    }

    // Process-specific: ABRICATE_SUMMARY
    withName: 'ABRICATE_SUMMARY' {
        ext.logs_subdir = "abricate-concat"
        ext.prefix = "abricate-${params.abricate_db}"
        ext.process_name = params.merge_folder
        ext.subdir = params.abricate_db
        ext.scope = "run"  // Override to run scope
    }
}
```

### Multi-Argument Configuration

For tools requiring separate argument sets (e.g., Kraken2 + Bracken):

```groovy
process {
    withName: 'BRACKEN' {
        // Primary arguments (for Kraken2)
        ext.args = [
            params.kraken2_quick_mode ? "--quick" : "",
            params.kraken2_use_mpa_style ? "--use-mpa-style" : "",
            "--confidence ${params.kraken2_confidence}",
            "--minimum-hit-groups ${params.kraken2_minimum_hit_groups}"
        ].join(' ').replaceAll("\\s{2,}", " ").trim()

        // Secondary arguments (for Bracken)
        ext.args2 = [
            "-l ${params.bracken_level}",
            "-t ${params.bracken_threshold}"
        ].join(' ').replaceAll("\\s{2,}", " ").trim()

        // Module-specific parameters
        ext.bracken_read_length = params.bracken_read_length
        ext.kraken2_keep_raw_output = params.kraken2_keep_raw_output
        ext.skip_krona = params.bracken_skip_krona
    }
}
```

## Usage in Modules

### Accessing Properties

```groovy
script:
// Access tool arguments
def args = task.ext.args ?: ''

// Access scope for output path construction
def scope = task.ext.scope

// Build conditional options
def extra_args = task.ext.args2 ?: ''
```

### Meta Construction Pattern

```groovy
script:
def prefix = task.ext.prefix ?: "${_meta.name}"
def meta = [:]
meta.id = "${prefix}-${task.process}"
meta.name = prefix
meta.scope = task.ext.scope
meta.process_name = task.ext.process_name

// Workflow-specific output paths
if (task.ext.wf == "teton") {
    meta.output_dir = "${prefix}/teton/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/teton/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
} else {
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
}
```

## Workflow Routing (task.ext.wf)

The `task.ext.wf` property controls output directory organization based on which workflow is running. This enables the same modules to be used by different entry workflows while maintaining distinct output structures.

### Workflow Values

| Value | Workflow | Output Structure |
|-------|----------|------------------|
| `"bactopia"` | Main Bactopia pipeline | `{sample}/tools/{process_name}/` |
| `"teton"` | Teton workflow | `{sample}/teton/tools/{process_name}/` |

### How Routing Works

1. **Configuration**: Set `ext.wf = params.wf` in module.config
2. **Module access**: Check `task.ext.wf` in the script section
3. **Path construction**: Conditionally build `meta.output_dir` based on workflow

### Output Path Examples

For a sample named "sample001" running the `mlst` tool:

**Bactopia workflow** (`task.ext.wf == "bactopia"`):

```bash
sample001/tools/mlst/sample001.tsv
sample001/tools/mlst/logs/mlst-run.log
```

**Teton workflow** (`task.ext.wf == "teton"`):

```bash
sample001/teton/tools/mlst/sample001.tsv
sample001/teton/tools/mlst/logs/mlst-run.log
```

### Implementation Pattern

In `module.config`:

```groovy
process {
    withName: 'MLST' {
        ext.wf = params.wf  // Inherits workflow from params
        // ... other config
    }
}
```

In `main.nf`:

```groovy
script:
def prefix = task.ext.prefix ?: "${_meta.name}"
def meta = [:]

// Construct output path based on workflow
if (task.ext.wf == "teton") {
    meta.output_dir = "${prefix}/teton/tools/${task.ext.process_name}/${task.ext.subdir}"
} else {
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
}
```

### Why This Pattern Exists

- **Separation of outputs**: Keeps Teton-specific results distinct from main Bactopia outputs
- **Module reuse**: Same module code works for both workflows
- **Configuration-driven**: Behavior controlled by config, not hardcoded
- **Clear provenance**: Output path indicates which workflow produced the results

### Scope-Based Output Paths

```groovy
// Sample scope (per-sample outputs)
if (task.ext.scope == "sample") {
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
}

// Run scope (aggregated outputs)
if (task.ext.scope == "run") {
    meta.output_dir = "${task.ext.process_name}"
    meta.logs_dir = "${task.ext.process_name}/logs/${task.ext.logs_subdir}/${task.ext.subdir}"
}
```

## Required vs Optional Properties

### Always Required

- `toolName` - For container/environment identification
- `docker` - Docker container reference
- `image` - Singularity image URL
- `condaDir` - Conda directory path
- `wf` - Workflow identifier
- `scope` - Output scope
- `process_name` - Process name for paths

### Usually Required

- `args` - Tool arguments (can be empty string)
- `subdir` - Output subdirectory (can be empty string)
- `logs_subdir` - Logs subdirectory (can be empty string)

### Optional

- `args2` - Secondary arguments (only for multi-step tools)
- `prefix` - Filename prefix override (defaults to `meta.name`)
- Module-specific properties (as needed)

## Naming Conventions

### toolName Format

Convert conda package specification to directory-safe name:

```groovy
ext.toolName = "bioconda::abricate=1.0.1".replace("=", "-").replace(":", "-").replace(" ", "-")
// Result: "bioconda-abricate-1.0.1"
```

### Argument Construction

Use array join pattern for clean, readable arguments:

```groovy
ext.args = [
    "--flag1",
    params.optional ? "--optional" : "",
    "--param ${params.value}"
].join(' ').replaceAll("\\s{2,}", " ").trim()
```

## See Also

- [Technical Specifications](../standards/03-technical-specs.md) - For type conventions
- [Module Documentation](../standards/05-module-documentation.md) - For module structure
- [Configuration System](../project/03-configuration-system.md) - For configuration hierarchy
