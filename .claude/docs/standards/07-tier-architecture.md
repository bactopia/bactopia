# Tier Architecture Standards

## Overview

Bactopia follows a three-tier architecture with strict separation of concerns. Each tier has defined responsibilities, allowed operations, and conventions. These rules are formalized to ensure consistency for maintenance and to enable reliable code generation by LLM agents.

## Tier 1: Workflows ("API Client")

Workflows are user-facing entry points. Users interact with workflows, and workflows are built exclusively from subworkflows.

### Rules
- Call ONLY subworkflows, never modules directly
- Channel manipulation uses named plugin functions from nf-bactopia, not raw operators
- Conditional routing (`if/else`) is allowed for selecting between subworkflows
- Cross-subworkflow synchronization (`.collect()`, `.map()`) is allowed when needed (e.g., collecting per-sample outputs for a run-level subworkflow)
- Publish/output blocks follow the standardized template
- Output aggregation with `.mix()` is allowed for combining subworkflow outputs

### Allowed Operations

| Operation | How | Example |
|-----------|-----|---------|
| Call subworkflows | Direct invocation | `ABRICATE(assembly)` |
| Aggregate outputs | `.mix()` | `ch_out = SUB1.out.mix(SUB2.out)` |
| Collect nf_logs | Plugin function | `collectNextflowLogs(ch_sample_outputs)` |
| Collect field for run-level input | Plugin function | `gather(ch_samples, 'gff', [name: 'panaroo'])` |
| Conditional routing | `if/else` | `if (!params.skip_phylogeny) { IQTREE(input) }` |
| Publish outputs | Standardized output block | See template below |

### Not Allowed in Workflows
- Calling modules directly
- Inline `.map()`, `.filter()`, `.flatMap()` for data transformation (use plugin functions)
- Channel operators for field extraction or record reshaping (belongs in subworkflows)

### Publish/Output Template

All workflows use this standardized structure:

```groovy
publish:
sample_outputs = ch_sample_outputs
sample_nf_logs = ch_sample_nf_logs
run_outputs    = ch_run_outputs       // optional
run_nf_logs    = ch_run_nf_logs       // optional

output {
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
```

### Workflow Types

| Type | Location | Has ext | Examples |
|------|----------|---------|---------|
| Named workflow | `workflows/<name>/` | No | bactopia, staphopia, teton, cleanyerreads |
| Bactopia tool | `workflows/bactopia-tools/<name>/` | Yes | abricate, pangenome, snippy |

Named workflows manage their own input gathering (via BACTOPIA_INIT). Bactopia tools use BACTOPIATOOL_INIT with `params.workflow.ext` to collect inputs from prior Bactopia runs.

---

## Tier 2: Subworkflows ("API Request")

Subworkflows handle all internal plumbing needed to connect modules and orchestrate analysis. This is where channel operators, field extraction, filtering, and aggregation logic belongs.

### Rules
- Every module has a corresponding subworkflow (1:1 wrapper rule)
- Handle all internal plumbing: field extraction, filtering, gather, concat
- Emit `sample_outputs` and/or `run_outputs` with standardized record structure
- Accept the most natural upstream output format
- May call other subworkflows (nesting is allowed)

### Allowed Operations

| Operation | How | Example |
|-----------|-----|---------|
| Call modules | Direct invocation | `ABRICATE_RUN(assembly)` |
| Call subworkflows | Direct invocation | `SNPDISTS(input)` |
| Aggregate per-sample results | `gatherCsvtk()` | `gatherCsvtk(MODULE.out, 'tsv', [name: 'tool'])` |
| Filter records | `filterWithData()` | `filterWithData(MODULE.out, ['fna', 'faa'])` |
| Channel operators | Any Nextflow operator | `.map()`, `.filter()`, `.combine()`, etc. |
| Conditional logic | `if/else` | Tool selection based on params |

### Output Convention

All subworkflows emit standardized outputs:

```groovy
emit:
sample_outputs: Channel<Record> = MODULE.out
run_outputs: Channel<Record>    = CSVTK_CONCAT.out
```

Record fields follow the module output convention: named fields for downstream use, plus standard generic fields (results, logs, nf_logs, versions) for publishing.

### Subworkflow Categories

| Category | Description | Examples |
|----------|-------------|---------|
| Sample wrapper | 1:1 module wrapper with optional gather/concat | abricate, amrfinderplus, busco |
| Run aggregator | Consumes collected samples, produces run-level output | pangenome, snippy/core, iqtree |
| Gateway | Routes to other subworkflows based on conditions | merlin |
| Infrastructure | Core pipeline mechanics | gather, qc, assembler, utils |

Gateway and infrastructure subworkflows are documented exceptions to the standard patterns.

---

## Tier 3: Modules ("API Server")

Modules are the most basic building blocks. This is where processing happens and data is generated. Modules have no awareness of upstream or downstream context.

### Rules
- Emit a single record with named downstream fields + standard generic fields
- No awareness of what calls them or what consumes their output
- Configuration via `ext` properties in module.config
- Do not perform channel manipulation

### Output Record Structure

Every module emits a record with two categories of fields:

```groovy
record(
    // Named fields -- used downstream by subworkflows
    meta: meta,
    tsv: file("${prefix}.tsv"),
    report: file("${prefix}-report.txt", optional: true),

    // Standard generic fields -- used for publishing
    results: [ files("${prefix}.tsv"), files("${prefix}-report.txt", optional: true) ],
    logs: files("*.{log,err}", optional: true),
    nf_logs: files(".command.*"),
    versions: files("versions.yml")
)
```

**Standard generic fields** (present in every module):
- `results` -- list of all publishable output files
- `logs` -- tool-specific log files
- `nf_logs` -- Nextflow command files
- `versions` -- software version tracking

**Named fields** are tool-specific and documented via GroovyDoc `@output`.

### Configuration via module.config

Each module has a `module.config` file that sets `ext` properties:

| Property | Purpose | Example |
|----------|---------|---------|
| `ext.wf` | Workflow identifier | `params.wf` |
| `ext.scope` | Output scope | `"sample"` or `"run"` |
| `ext.process_name` | Process name for output paths | `"abricate"` |
| `ext.subdir` | Output subdirectory | `""` or `params.run_name` |
| `ext.args` | Tool CLI arguments | Assembled from params |
| `ext.toolName` | Conda package spec | `"bioconda::abricate=1.0.1"` |
| `ext.docker` | Docker image URL | `"biocontainers/abricate:1.0.1--..."` |
| `ext.image` | Singularity image URL | `"https://depot.galaxyproject.org/..."` |

---

## Plugin Functions

Channel manipulation is standardized through named functions in the nf-bactopia plugin. This keeps workflows and subworkflows readable and centralizes operator logic.

### Function Reference

| Function | Purpose | Used in |
|----------|---------|---------|
| `gather(ch, field, meta)` | Extract field, collect across samples, return record | Workflows (run-level input prep) |
| `gatherCsvtk(ch, field, meta)` | Same as gather, output field renamed to `csv` | Subworkflows (CSVTK_CONCAT input) |
| `gatherFields(ch, fieldMapping, meta)` | Multi-field gather with rename map | Future use |
| `filterWithData(ch, fields)` | Filter records by non-null fields, project down | Subworkflows |
| `collectNextflowLogs(ch)` | Extract nf_logs from records into [meta, file] tuples | All workflows |

All gather functions delegate to a private `_gather(ch, Map fieldMapping, Map meta)` core.

---

## Extension System

Bactopia tools declare their input requirements via `params.workflow.ext` in their `nextflow.config`.

### Format

`params.workflow.ext` is a list of keys from the controlled vocabulary:

```groovy
params {
    workflow {
        name = "amrfinderplus"
        ext = ["fna", "faa", "gff"]
    }
}
```

### Controlled Vocabulary

| Key | Meaning |
|-----|---------|
| `fna` | Assembled genome |
| `fna_anno` | Annotator-formatted assembly |
| `faa` | Protein sequences |
| `gff` | Gene coordinates |
| `gbk` | GenBank file |
| `tsv_meta` | Analysis metadata TSV |
| `blastdb` | BLAST database |
| `r1` | Forward reads (Illumina PE) |
| `r2` | Reverse reads (Illumina PE) |
| `se` | Single-end reads (Illumina SE) |
| `lr` | Long reads (ONT/PacBio) |
| `fastq` | Alias: expands to r1, r2, se, lr |

### Semantics
- Non-read keys are **all required** (AND logic): sample must have every listed file
- Read keys declare **accepted read types**: sample needs at least one valid read set
- A sample will not have both PE (r1+r2) and SE -- they are mutually exclusive
- A sample CAN have both short reads and LR (hybrid)
- `fastq` is an alias that expands to all four read keys

### Examples

```groovy
["fna"]                    // Assembly only (abricate, mlst)
["fna", "faa", "gff"]     // Assembly + annotation (amrfinderplus)
["fna", "tsv_meta"]       // Assembly + metadata (quast)
["r1", "r2", "se"]        // Illumina reads, PE or SE (snippy)
["r1", "r2"]              // Illumina PE only (seroba, pneumocat)
["fastq"]                 // Any reads (ariba, kraken2)
["fna", "fastq"]          // Assembly + any reads (merlin)
["faa"]                   // Proteins only (eggnog)
["blastdb"]               // BLAST database (blastn, blastp)
["gbk"]                   // GenBank file (phispy)
["gff"]                   // GFF file (pangenome)
```

### Validation
- A lint rule enforces that every bactopia-tool nextflow.config has `params.workflow.ext` as a list
- Every value must be from the controlled vocabulary
- Named workflows (bactopia, staphopia, teton, cleanyerreads) do not use ext

---

## catalog.json

A machine-readable index of all workflows, subworkflows, and modules. Auto-generated by the `bactopia-catalog` CLI command from bactopia-py. Replaces `data/workflows.yml`.

### Data Sources

All data is derived from the codebase -- nothing is hand-maintained:
- Descriptions from GroovyDoc
- Dependencies from `include` statements
- Paths from directory structure
- Workflow type from path convention
- ext values from nextflow.config
- Takes/emits contracts from record types
- Tool versions from `ext.toolName` in module.config
- Scope from `ext.scope` in module.config

### Usage
- LLM agents load catalog.json to understand component contracts and composition
- bactopia-py CLI commands use it for environment downloads, schema merging, etc.
- Regenerate after adding or modifying components: `bactopia-catalog`
