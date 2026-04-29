# Technical Specifications

## Overview
This document contains the technical specifications, conventions, and "gotchas" for implementing components in the Bactopia pipeline.

## Static Typing Conventions

### Universal Requirements
- ALL Nextflow files must have `nextflow.enable.types = true` at the top (required while typed records are in preview)
- Static typing is fully implemented across all components
- No files should be missing type annotations

### File Output Types

#### Single Files: `file()`
- **Returns**: `Path`
- **Used in**: Named record fields for downstream access
- **Example**: `tsv: file("${prefix}.tsv")` (inside `record()` output block)
- **Use cases**: Single, known files with predictable names

#### Multiple Files: `files()`
- **Returns**: `Set<Path>`
- **Used in**: Generic record fields and `results` lists
- **Example**: `logs: files("*.{log,err}", optional: true)` (inside `record()` output block)
- **Use cases**: Wildcards, optional files, or multiple files

**Key Rule**: Use `file()` for single, known files. Use `files()` for wildcards, optional files, or multiple files.

### Channel Type Declarations

Subworkflows declare channel types on `take:` block entries as `Channel<Record>`:

```nextflow
workflow ABRICATE {
    take:
    assembly: Channel<Record>

    main:
    ch_abricate_run = ABRICATE_RUN(assembly)
    ...
}
```

Intermediate channels (e.g. `ch_abricate_run` above) inherit their `Channel<Record>` type from the module or subworkflow they are assigned from — no explicit `as Channel<...>` cast is needed.

### Standard Module Input Types

The first positional input is a typed `record (...)` block. Additional inputs are declared on separate lines.

**Single assembly** (e.g. [modules/prokka/main.nf:45-51](../../../modules/prokka/main.nf#L45-L51)):

```nextflow
input:
record (
    meta: Record,
    fna: Path
)
proteins   : Path?
prodigal_tf: Path?
```

**Multi-read inputs** — modules accepting any mix of read types (e.g. [modules/bactopia/assembler/main.nf:51-59](../../../modules/bactopia/assembler/main.nf#L51-L59)):

```nextflow
input:
record (
    meta: Record,
    r1: Path?,
    r2: Path?,
    se: Path?,
    lr: Path?,
    fna: Path?
)
```

- `r1`: Illumina paired-end forward
- `r2`: Illumina paired-end reverse
- `se`: Single-end Illumina reads
- `lr`: Long reads (ONT/PacBio)

**Multiple distinct inputs** — reads plus a reference (e.g. [modules/snippy/run/main.nf:58-65](../../../modules/snippy/run/main.nf#L58-L65)):

```nextflow
input:
record (
    meta: Record,
    r1: Path?,
    r2: Path?,
    se: Path?
)
reference: Path
```

**Note**: Input blocks use named `Record` fields rather than positional `Tuple<>` annotations. In the module `script:` block, `def _meta = meta` aliases the input record before `meta = record(...)` constructs a new output record.

## Path? Optional Parameters

### Declaration

Optional file parameters use Nextflow's native `Path?` type. Mark the field with `?` inside the `record (...)` block when the slot itself is optional, or use `Path?` for additional inputs declared on their own line:

```groovy
input:
record (
    meta: Record,
    fna: Path
)
proteins   : Path?
prodigal_tf: Path?
```

When a `Path?` parameter has no value, it is `null`.

### Detection in Scripts

Use null checks to conditionally build command-line arguments:

```groovy
def proteins_opt = proteins != null ? "--proteins ${proteins.getName()}" : ""
```

### GroovyDoc Convention

Optional fields are marked with a `?` suffix in GroovyDoc to mirror the `Path?` type:

```groovy
 * @input record(meta, r1?, r2?, se?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 *
 * @input proteins?
 * FASTA file of trusted proteins to first annotate from
```

Rules:
- Add `?` to any field backed by `Path?` in the take block
- Add `?` to output record fields that use `optional: true` in `file()`
- **Never** add `?` to standard fields: `meta`, `results`, `logs`, `nf_logs`, `versions`

### Output with Optional Fields

```groovy
output:
record(
    meta: meta,
    r1: r1 != null ? file("fastqs/${prefix}_R1.fastq.gz", optional: true) : null,
    fna: file("assembly/${prefix}.fna.gz"),
)
```

Fields using `optional: true` or conditional null passthrough get `?` in the `@output record(...)` line.

## Channel Output Patterns

### Module Output Pattern (record block)

Modules emit a single `record()` with named fields (for downstream access) and generic fields (for publishing):

```nextflow
output:
record(
    // Named fields (used downstream)
    meta: meta,
    tsv: file("${prefix}.tsv"),
    // Generic fields (used for publishing)
    results: [
        files("${prefix}.tsv")
    ],
    logs: files("*.{log,err}", optional: true),
    nf_logs: files(".command.*"),
    versions: files("versions.yml")
)
```

- Use `file()` for named fields that will be accessed downstream (returns `Path`)
- Use `files()` in the `results` list and for logs/versions (returns `Set<Path>`)
- Use exactly **one space** after the colon in record fields (e.g., `meta: meta`, not `meta:  meta`)

### Subworkflow Output Pattern (2 emit channels)

```nextflow
workflow ABRICATE {
    take:
    assembly: Channel<Record>

    main:
    ch_abricate_run = ABRICATE_RUN(assembly)
    ch_abricate_summary = ABRICATE_SUMMARY(gatherFields(ch_abricate_run, [tsv: 'reports'], [name: 'abricate']))

    emit:
    sample_outputs = ch_abricate_run
    run_outputs = ch_abricate_summary
}
```

Module results are assigned directly to intermediate variables (`ch_module = MODULE(...)`) and those variables are emitted — the old `MODULE.out` form is no longer used. Subworkflows emit two channels:

- `sample_outputs`: The per-sample module record (passed through directly)
- `run_outputs`: Aggregated results from a summary module (typically `CSVTK_CONCAT` or similar)

## Meta Record Structure

The meta record is a `Record` that carries sample metadata through the pipeline. It is passed as the first field of input/output records and contains both input-derived and runtime-constructed properties.

### Complete Schema

#### Input-Derived Properties

These properties come from sample sheets or input initialization:

| Property | Type | Description | Example |
|----------|------|-------------|---------|
| `id` | String | Unique sample identifier | `"sample001"` |
| `name` | String | Sample name (often same as id) | `"sample001"` |
| `runtype` | String | Read type classification | `"paired-end"`, `"single-end"`, `"ont"`, `"hybrid"` |
| `single_end` | Boolean | True if single-end reads | `true`, `false` |
| `genome_size` | Integer | Estimated genome size in bp | `5000000` |
| `species` | String | Species name for analysis | `"Staphylococcus aureus"` |

#### Runtime-Constructed Properties

These properties are set within module scripts based on `task.ext` configuration:

| Property | Type | Description | Set From |
|----------|------|-------------|----------|
| `id` | String | Process-specific identifier | `"${prefix}-${task.process}"` |
| `name` | String | Sample name prefix | `task.ext.prefix ?: _meta.name` |
| `scope` | String | Output scope | `task.ext.scope` |
| `output_dir` | String | Output directory path | Constructed from task.ext |
| `logs_dir` | String | Logs directory path | Constructed from task.ext |
| `process_name` | String | Current process name | `task.ext.process_name` |

#### Module-Specific Properties

Some modules add additional properties:

| Property | Type | Modules | Description |
|----------|------|---------|-------------|
| `is_compressed` | Boolean | qc, gather | Input files are gzipped |
| `original_runtype` | String | qc, gather | Runtype before normalization |
| `teton_reads` | Boolean | teton workflow | Using Teton read processing |

### Standard Meta Construction Pattern

```groovy
script:
def _meta = meta
prefix = task.ext.prefix ?: "${_meta.name}"

// Create a new meta record for this process
meta = record(
    id: "${prefix}-${task.process}",
    name: prefix,
    scope: task.ext.scope,
    output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
    logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
    process_name: task.ext.process_name
)
```

Modules that need to emit under a different workflow path (e.g. `teton`) introduce a `wfPath` helper and interpolate it into the `output_dir`/`logs_dir` fields:

```groovy
def String wfPath = task.ext.wf == "teton" ? "teton/main" : "main"

meta = record(
    // ... other fields ...
    output_dir: "${prefix}/${wfPath}/${task.ext.process_name}/${task.ext.subdir}",
    logs_dir: "${prefix}/${wfPath}/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
    process_name: task.ext.process_name
)
```

### Output Directory Patterns

Records are immutable, so `output_dir` and `logs_dir` must be set inside the `meta = record(...)` constructor shown above — never as post-construction assignments. The two standard forms are:

```groovy
// Sample scope (inside meta = record(...))
output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
logs_dir:   "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",

// Run scope (inside meta = record(...))
output_dir: "${task.ext.process_name}",
logs_dir:   "${task.ext.process_name}/logs/${task.ext.logs_subdir}/${task.ext.subdir}",
```

### Scope Values

| Scope | Description | Output Pattern |
|-------|-------------|----------------|
| `sample` | Per-sample outputs | `{sample}/tools/{process}/{subdir}/` |
| `run` | Aggregated outputs | `{process}/` |

### Runtype Values

| Runtype | Description |
|---------|-------------|
| `paired-end` | Illumina paired-end reads (R1 + R2) |
| `single-end` | Illumina single-end reads |
| `ont` | Oxford Nanopore long reads |
| `hybrid` | Illumina + long reads combined |
| `merge-pe` | Merged paired-end reads |
| `merge-se` | Merged single-end reads |

## Variable Naming Conventions

### Two-Level Naming Convention

Record field names and channel names follow a dual convention:
- **Record fields** = bioinformatics format abbreviation (short, technical)
- **Channel names** (subworkflow takes/emits, workflow variables) = human-readable concept

**Singular vs plural**:
- Singular = per-sample (one record per sample)
- Plural = aggregated / collection intended for multi-sample analysis

#### Canonical Name Table

| Concept | Record field | Channel (singular) | Channel (plural) |
|---|---|---|---|
| Assembly | `fna` | `assembly` | `assemblies` |
| Proteins | `faa` | `proteins` | `proteins` |
| Feature nt FASTA | `ffn` | (rarely channeled) | |
| GFF annotation | `gff` | `gff` | `gffs` |
| GenBank annotation | `gbff` | `gbff` | `gbffs` |
| MSA / alignment | `aln` | `alignment` | `alignments` |
| Filtered alignment | `filtered_aln` | | |
| Masked alignment | `masked_aln` | | |
| Full alignment | `full_aln` | | |
| Clean full alignment | `clean_full_aln` | | |
| Phylogenetic tree | `nwk` | `tree` | `trees` |
| Per-sample ref-aligned FA | `aligned_fa` | (keep as-is) | |
| Variant calls | `vcf` | `vcf` | `vcfs` |
| Distance matrix | `dist` | `dist` | |
| Results (tab) | `tsv` | `tsv` | |
| Results (comma) | `csv` | `csv` | |

#### Examples

Module output uses format abbreviation for record fields:
```nextflow
record(meta: meta, aln: file("${prefix}.aln.gz"), nwk: file("${prefix}.treefile"), ...)
```

Subworkflow takes use human-readable channel names:
```nextflow
take:
alignment: Channel<Record>  // contains records with field `aln`
```

Subworkflow emits use human-readable names:
```nextflow
emit:
alignment = MODULE.out.map { r -> record(meta: ..., aln: r.masked_aln) }
```

GroovyDoc `@input`/`@output` documents the record field names (format abbreviation):
```
@input record(meta, aln)
@output record(meta, aln, nwk, results, logs, nf_logs, versions)
```

### Channel Naming
- Use descriptive names: `ch_results`, `ch_logs`, `ch_nf_logs`, `ch_versions`
- For component-specific outputs: use the output type (e.g., `aln`, `csv`, `report`)

## Explicit Positional Read Tuple Pattern

Modules that process sequencing reads use a 5-slot positional tuple pattern to handle different read types explicitly.

### The Pattern

```nextflow
// Input signature
input:
record (
    meta: Record,
    r1: Path?,
    r2: Path?,
    se: Path?,
    lr: Path?
)
```

| Slot | Variable | Description |
|------|----------|-------------|
| 1 | `meta` | Sample metadata record |
| 2 | `r1` | Illumina R1 (paired-end forward) |
| 3 | `r2` | Illumina R2 (paired-end reverse) |
| 4 | `se` | Single-end Illumina reads |
| 5 | `lr` | Long reads (ONT/PacBio) |

### Pre-GATHER vs Post-GATHER

The GATHER module transforms read channels. In both stages the channel element is a `Record` containing a `meta: Record` and per-slot read fields — what changes is how many files each slot holds:

| Stage | Channel element | Purpose |
|-------|-----------------|---------|
| Pre-GATHER | `Channel<Record>` with `record(meta, r1_files: Set<Path?>, r2_files: Set<Path?>, se_files: Set<Path?>, lr_files: Set<Path?>, fna_files: Set<Path?>)` | Multiple files per slot (multi-lane/run, pre-merge) |
| Post-GATHER | `Channel<Record>` with `record(meta, r1?, r2?, se?, lr?, fna?)` where each read slot is `Path?` | Single consolidated file per slot (post-merge) |

### Detecting Read Type in Scripts

```groovy
script:
// Check which slots are populated
has_r1 = r1 != null
has_r2 = r2 != null
has_se = se != null
has_lr = lr != null

// Determine read type
meta.single_end = has_se && !has_r1 && !has_r2
meta.is_paired = has_r1 && has_r2

// Build tool-specific read arguments
def read_inputs = has_lr ? "${lr}" : (meta.single_end ? "${se}" : "${r1} ${r2}")
```

### Modules Using This Pattern

- `bracken`, `kraken2` - Taxonomic classification
- `snippy/run` - Variant calling
- `tbprofiler/profile` - TB profiling
- `mykrobe/predict` - AMR prediction
- `ariba/run` - ARIBA analysis (paired-end only)
- `bactopia/qc` - Quality control
- `bactopia/gather` - Read consolidation

## stageAs Directive

The `stageAs` directive organizes input files into specific directory structures within a process's working directory.

### Purpose

When inputs are `Set<Path>` collections or need organized staging, `stageAs` creates predictable directory structures that scripts can reference reliably.

### Basic Pattern

```nextflow
stage:
    stageAs <variable_name>, '<directory_pattern>'
```

The variable comes **first**, then the staging path. Current modules use a `staging/<input_type>/*` convention, where `<input_type>` matches the record field name (`fna`, `r1`, `r2`, `se`, `lr`, `gff`, `faa`, `csv`, `json`, etc.). This keeps scripts predictable — `staging/fna/*` always means assemblies, `staging/r1/*` always means R1 reads.

### Common Use Cases

#### Single Collection

```nextflow
// modules/roary/main.nf - pangenome GFF files
stage:
    stageAs gff, 'staging/gff/*'
```

#### Multiple Input Types

```nextflow
// modules/stecfinder/main.nf - mixed reads + assembly
stage:
    stageAs fna, 'staging/fna/*'
    stageAs r1,  'staging/r1/*'
    stageAs r2,  'staging/r2/*'
    stageAs se,  'staging/se/*'
    stageAs lr,  'staging/lr/*'
```

#### Assembly + Database

```nextflow
// modules/gtdbtk/classifywf/main.nf
stage:
    stageAs fna, 'staging/fna/*'
    stageAs db,  'staging/gtdb/*'
```

#### Multi-Lane Reads (GATHER only)

```nextflow
// modules/bactopia/gather/main.nf - pre-merge, multiple files per slot
stage:
    stageAs r1_files,  "staging/r1/*???-r1"
    stageAs r2_files,  "staging/r2/*???-r2"
    stageAs se_files,  "staging/se/*???-se"
    stageAs lr_files,  "staging/lr/*???-lr"
    stageAs fna_files, "staging/fna/*???-assembly"
```

GATHER is the only module that appends a `*???-<type>` suffix — it runs before per-slot merging, so the suffix disambiguates multiple files in a single slot (e.g. lane-split R1 files). Post-GATHER modules always have one file per slot and use the plain `staging/<type>/*` form.

### Modules Using stageAs

| Module(s) | Pattern | Purpose |
|-----------|---------|---------|
| `bactopia/gather` | `staging/{r1,r2,se,lr,fna}/*???-{type}` | Multi-lane read consolidation (pre-merge) |
| `bactopia/qc`, `stecfinder`, `amrfinderplus/run` | `staging/{r1,r2,se,lr,fna,...}/*` | Mixed reads + assembly inputs |
| `roary`, `pirate`, `panaroo/run` | `staging/gff/*` | GFF collection for pangenome tools |
| `csvtk/concat` | `staging/csv/*` | TSV/CSV file collection for aggregation |
| `gtdbtk/classifywf` | `staging/fna/*`, `staging/gtdb/*` | Assembly + database staging |
| `tbprofiler/collate` | `staging/json/*` | JSON result collection |
| `prokka`, `agrvate`, `bakta/run` | `staging/fna/*` | Single assembly staging |
| `defensefinder/run` | `staging/faa/*` | Protein FASTA staging |

## Database Handling Patterns

Many modules accept external databases. The codebase supports two formats:
- **Directory**: For local/HPC systems (faster access)
- **Tarball (.tar.gz)**: For cloud systems (easier transport)

### Tarball Detection Pattern

```groovy
script:
def is_tarball = db.getName().endsWith(".tar.gz") ? true : false

// In shell:
if [ "${is_tarball}" == "true" ]; then
    mkdir database
    tar -xzf ${db} -C database
    DB_PATH=$(find database/ -name "hash.k2d" | sed 's=hash.k2d==')
else
    DB_PATH="${db}"
fi
```

### Database Path Discovery

Different databases have different marker files:

| Database | Marker File | Discovery Command |
|----------|-------------|-------------------|
| Kraken2/Bracken | `hash.k2d` | `find database/ -name "hash.k2d"` |
| MLST | `mlst.db` | `find database/ -name "mlst.db"` |
| GTDB-Tk | `metadata.txt` | `find database/ -name "metadata.txt"` |
| Bakta | `bakta.db` | `find database/ -name "bakta.db"` |

### Modules with Database Support

| Module | Database Type | Notes |
|--------|---------------|-------|
| `bracken`, `kraken2` | Kraken2 DB | Tarball or directory |
| `mlst` | PubMLST | Tarball or directory |
| `gtdbtk/classifywf` | GTDB | Large database, tarball recommended for cloud |
| `bakta/run` | Bakta DB | Supports light/full versions |
| `checkm2/predict` | CheckM2 | Diamond database |
| `midas/species` | MIDAS | Species database |
| `blast/*` | BLAST DB | Always tarball extraction |

### Best Practice

```groovy
// Module config - allow both formats
params {
    tool_db = null  // User provides path to directory or .tar.gz
}

// In module script - handle both
def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
```

## Utility Functions

The `nf-bactopia` plugin exports a set of channel-manipulation helpers used throughout subworkflows. This section summarizes the ones most commonly needed when writing a new module or subworkflow — see [Plugin Functions](../reference/04-plugin-functions.md) for the full reference.

| Function | Purpose |
|----------|---------|
| `gather` | Collect a single named field from records into a `Set`, emit one aggregated record. |
| `gatherCsvtk` | Like `gather`, but renames the collected field to `csv` for direct feeding into `CSVTK_CONCAT`. |
| `gatherFields` | Collect multiple fields with explicit input→output field name mapping. |
| `filterWithData` | Drop records where every listed field is `null`; projects to `meta` + requested fields. |
| `combineWith` | Cartesian product between a single-item channel and a multi-item channel, merging each item under a named field. |
| `collectNextflowLogs` | Expand a record's `nf_logs` field into individual `[meta, file]` tuples for publishing. |

Typical usage in a subworkflow:

```nextflow
include { gatherFields } from 'plugin/nf-bactopia'

ch_abricate_run = ABRICATE_RUN(assembly)
ch_abricate_summary = ABRICATE_SUMMARY(gatherFields(ch_abricate_run, [tsv: 'reports'], [name: 'abricate']))
```

## Container Directives

All modules use the short-form container directive:

```groovy
conda "${task.ext.condaDir}/${task.ext.toolName}"
container "${task.ext.container}"
```

Container resolution (choosing Docker vs Singularity image) is handled globally in `conf/base.config` via a `withName: '.*'` block that sets `ext.container` based on the container engine. Do **not** use the inline ternary pattern in individual modules.

## Script Block Structure

Module script blocks must include a `# Cleanup` comment before the version tracking block. This marks where cleanup operations (e.g., removing temporary files) should go:

```bash
    # Cleanup
    rm -rf temp_dir/

    cat <<-END_VERSIONS > versions.yml
    ...
```

If no cleanup is needed, the comment is still required as a placeholder for consistency:

```bash
    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    ...
```

## Version Tracking
Always include a `versions.yml` file with software version information:

```nextflow
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    ${task.ext.toolName}: \$(echo \$(MODULE --version 2>&1) | sed 's/^.*MODULE //' )
END_VERSIONS
```

Expected format:

```yaml
"PROCESS_NAME":
    TOOL_NAME: VERSION_STRING
```

## Common Gotchas

### What NOT to Do
1. **Don't "fix" Path? workarounds** - They're necessary until Nextflow improves
2. **Don't change `file()` to `files()`** - The difference is intentional
3. **Don't modify typing without understanding** - Types are deliberately chosen
4. **Don't alter the subworkflow 2-channel emit pattern** (`sample_outputs`, `run_outputs`) - It's fundamental to the architecture

### Type Checking
- Validate with `nextflow config` to check type annotations
- Ensure consistent typing across connected components
- Remember: `file()` returns `Path`; `files()` returns `Set<Path>`. Choose based on whether the field holds one known file or a wildcard / optional / multi-file collection.

## See Also

- [Style Guide](../standards/01-style-guide.md) - For template formats
- [Logic Rules](../standards/02-logic-rules.md) - For classification logic
- [Plugin Functions](../reference/04-plugin-functions.md) - For gather() and related channel utilities
- [task.ext Properties](../reference/05-task-ext-properties.md) - For module configuration
- [Troubleshooting](../reference/02-troubleshooting.md) - For common error solutions
