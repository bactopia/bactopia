# Technical Specifications

## Overview
This document contains the technical specifications, conventions, and "gotchas" for implementing components in the Bactopia pipeline.

## Static Typing Conventions

### Universal Requirements
- ALL Nextflow files must have `nextflow.preview.types = true` at the top (until Nextflow v26.04)
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
Standard pattern for workflow-level channels:

```nextflow
ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>
```

### Standard Module Input Types

The codebase uses Record-typed inputs with explicit named parameters:

- **Single assemblies**: `(_meta: Map, assembly: Path): Record` - For modules processing a single assembly file
- **Multi-read inputs**: `(_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record` - For modules accepting multiple read types:
    - `r1`: Illumina paired-end forward
    - `r2`: Illumina paired-end reverse
    - `se`: Single-end Illumina reads
    - `lr`: Long reads (ONT/PacBio)
- **Multiple distinct inputs**: `(_meta: Map, assembly: Path, meta_file: Path): Record` - For modules requiring multiple files (e.g., assembly + metadata)
- **Additional inputs** are declared on separate lines: `db: Path`, `proteins: Path?`, etc.

**Note**: The codebase uses `Record` return types with named parameters (e.g., `_meta: Map`) rather than `Tuple<>` type annotations.

## Path? Optional Parameters (Workarounds)

### The Problem
Nextflow has incomplete support for optional Path parameters (`Path?`). This is a known bug being worked on by the Nextflow team.

### The Workaround: EMPTY_* Files
Bactopia uses placeholder files to indicate the absence of optional files:

```groovy
// Default parameter pattern
bakta_db : Path? = "${projectDir}/data/empty/EMPTY_DB"
```

### Detection in Scripts

```groovy
def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""
```

**Note**: Use `.getName()` directly on `Path?` parameters. The older `.toList()[0].getName()` pattern is deprecated.

### Available EMPTY_* Files
Located in `/data/empty/`:
- `EMPTY_ADAPTERS` - For adapter files
- `EMPTY_ASSEMBLY` - For assembly files
- `EMPTY_BLASTDB` - For BLAST database files
- `EMPTY_DB` - For database files
- `EMPTY_EXTRA` - For extra files
- `EMPTY_GBK` - For GenBank files
- `EMPTY_META` - For metadata files
- `EMPTY_ONT` - For Oxford Nanopore long reads
- `EMPTY_PHIX` - For PhiX files
- `EMPTY_PRODIGAL_TF` - For Prodigal training files
- `EMPTY_PROTEINS` - For protein files
- `EMPTY_R1` - For read 1 files (paired-end forward)
- `EMPTY_R2` - For read 2 files (paired-end reverse)
- `EMPTY_REPLICONS` - For replicon files
- `EMPTY_SE` - For single-end Illumina reads
- `EMPTY_TF` - For training files

**Important**: These are TEMPORARY workarounds. Do NOT attempt to "fix" them. They will be removed when Nextflow properly supports Path?.

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
emit:
sample_outputs = MODULE.out
run_outputs = CSVTK_CONCAT.out
```

Subworkflows emit two channels:
- `sample_outputs`: The module's record output (passed through directly)
- `run_outputs`: Aggregated results from CSVTK_CONCAT (or similar aggregation module)

## Meta Map Structure

The meta map is a Map object that carries sample metadata through the pipeline. It is passed as the first element of input/output tuples and contains both input-derived and runtime-constructed properties.

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

### Output Directory Patterns

```groovy
// Sample scope
meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"

// Run scope
meta.output_dir = "${task.ext.process_name}"
meta.logs_dir = "${task.ext.process_name}/logs/${task.ext.logs_subdir}/${task.ext.subdir}"
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

### Standard Renames
- `fasta` → **`assembly`**: "Assembled contigs in FASTA format"
- `fastq`/`reads` → **`reads`**: "FASTQ reads (Illumina or Nanopore)"

### Channel Naming
- Use descriptive names: `ch_results`, `ch_logs`, `ch_nf_logs`, `ch_versions`
- For component-specific outputs: use the output type (e.g., `aln`, `csv`, `report`)

## Explicit Positional Read Tuple Pattern

Modules that process sequencing reads use a 5-slot positional tuple pattern to handle different read types explicitly.

### The Pattern

```groovy
// Input signature
input:
(_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record
```

| Slot | Variable | Description |
|------|----------|-------------|
| 1 | `meta` | Sample metadata map |
| 2 | `r1` | Illumina R1 (paired-end forward) |
| 3 | `r2` | Illumina R2 (paired-end reverse) |
| 4 | `se` | Single-end Illumina reads |
| 5 | `lr` | Long reads (ONT/PacBio) |

### Pre-GATHER vs Post-GATHER

The GATHER module transforms read channels:

| Stage | Type | Purpose |
|-------|------|---------|
| Pre-GATHER | `Tuple<Map, Set<Path>, Set<Path>, Set<Path>, Set<Path>>` | Multiple files per slot (multi-lane/run) |
| Post-GATHER | `Tuple<Map, Path?, Path?, Path?, Path?>` | Single consolidated file per slot |

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

```groovy
stage:
    stageAs '<directory_pattern>', <variable_name>
```

### Common Use Cases

#### Organizing Reads by Type

```groovy
// In modules/bactopia/gather/main.nf
stage:
    stageAs '*???-r1', r1_files
    stageAs '*???-r2', r2_files
    stageAs '*???-se', se_files
    stageAs '*???-lr', lr_files
```

The `*???-r1` pattern creates unique names that `find -name "*-r1"` can locate.

#### Staging File Collections

```groovy
// In modules/roary/main.nf - for pangenome GFF files
stage:
    stageAs 'gff-tmp/*', gff

// In script:
mkdir gff
cp -L gff-tmp/* gff/
```

#### Organizing Multiple Input Types

```groovy
// In modules/stecfinder/main.nf
stage:
    stageAs 'fna/*', fna
    stageAs 'reads/r1/*', r1
    stageAs 'reads/r2/*', r2
    stageAs 'reads/se/*', se
    stageAs 'reads/lr/*', lr
```

### Modules Using stageAs

| Module | Pattern | Purpose |
|--------|---------|---------|
| `bactopia/gather` | `*???-r1`, etc. | Multi-lane read organization |
| `roary`, `pirate`, `panaroo/run` | `gff-tmp/*` | GFF collection staging |
| `csvtk/concat` | `inputs/*` | CSV file collection |
| `gtdbtk/classifywf` | `fna-tmp/*`, `gtdb/*` | Assembly + database staging |
| `tbprofiler/collate` | `results-tmp/*` | JSON result collection |
| `prokka`, `agrvate` | `input/*` | Single assembly staging |

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

### gather
Collects a specific field from module record outputs and merges them for aggregation:

```nextflow
include { gather } from 'plugin/nf-bactopia'
CSVTK_CONCAT(gather(MODULE.out, 'tool-name', field: 'tsv'), 'tsv', 'tsv')
```

The `field:` named parameter extracts the specified field from each record in the channel.

### flattenPaths (deprecated)

`flattenPaths` was previously used in subworkflows to convert `Tuple<Map, Set<Path>>` to `Tuple<Map, Path>`. It is no longer used in subworkflows — they now pass through module record outputs directly via `sample_outputs` and `run_outputs`.

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
- Remember: `Tuple<Map, Set<Path>>` vs `Tuple<Map, Path>` is based on `files()` vs `file()`

## See Also

- [Style Guide](../standards/01-style-guide.md) - For template formats
- [Logic Rules](../standards/02-logic-rules.md) - For classification logic
- [Plugin Functions](../reference/04-plugin-functions.md) - For gather() and flattenPaths() usage
- [task.ext Properties](../reference/05-task-ext-properties.md) - For module configuration
- [Troubleshooting](../reference/02-troubleshooting.md) - For common error solutions
