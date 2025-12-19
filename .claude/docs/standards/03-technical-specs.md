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
- **Used in**: `Tuple<Map, Path>`
- **Example**: `versions = tuple(meta, file("versions.yml"))`
- **Use cases**: Single, known files with predictable names

#### Multiple Files: `files()`
- **Returns**: `Set<Path>`
- **Used in**: `Tuple<Map, Set<Path>>`
- **Example**: `logs = tuple(meta, files("*.{log,err}", optional: true))`
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

The codebase uses explicit positional types for clarity:

- **Single assemblies**: `Tuple<Map, Path>` - For modules processing a single assembly file
- **Multi-read inputs**: `Tuple<Map, Path?, Path?, Path?, Path?>` - For modules accepting multiple read types:
    - Position 1: R1 (Illumina paired-end forward)
    - Position 2: R2 (Illumina paired-end reverse)
    - Position 3: SE (Single-end Illumina)
    - Position 4: LR (Long reads - ONT/PacBio)
- **Multiple distinct inputs**: `Tuple<Map, Path, Path>` - For modules requiring multiple files (e.g., assembly + metadata)

**Note**: The codebase has transitioned from `Set<Path>` to explicit `Path` types for module inputs to provide clearer type safety and documentation.

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
- `EMPTY_DB` - For database files
- `EMPTY_EXTRA` - For extra files
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

### Module Output Pattern (3 channels)

```nextflow
output:
logs        = tuple(meta, files("*.{log,err}", optional: true))
nf_logs     = tuple(meta, files(".command.*"))
versions    = tuple(meta, files("versions.yml"))
```

**Note**: All standard outputs use `files()` for consistency, even for single files like `versions.yml`.

Modules may emit additional output channels as needed.

### Subworkflow Output Pattern (4 channels)

```nextflow
emit:
results    = flattenPaths([ch_results])
logs       = flattenPaths([ch_logs])
nf_logs    = flattenPaths([ch_nf_logs])
versions   = flattenPaths([ch_versions])
```

Subworkflows must always emit these four standard channels.

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
(_meta, r1, r2, se, lr) : Tuple<Map, Path?, Path?, Path?, Path?>
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

### flattenPaths
Converts `Tuple<Map, Set<Path>>` to `Tuple<Map, Path>` for outputs:

```nextflow
include { flattenPaths } from 'plugin/nf-bactopia'
results = flattenPaths([ch_results])
```

### gather
Merges multiple channels for aggregate operations:

```nextflow
include { gather } from 'plugin/nf-bactopia'
SUMMARY(gather(INPUT.out.report, 'tool-name'))
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
4. **Don't alter the subworkflow 4-channel pattern** - It's fundamental to the architecture

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
