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

Standard meta fields used across all components:

```groovy
meta.id = "${prefix}-${task.process}"
meta.name = prefix
meta.scope = task.ext.scope          // 'sample' or 'run'
meta.output_dir = "..."               // Used in output blocks
meta.logs_dir = "..."                 // Used in output blocks
meta.process_name = task.ext.process_name
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

## Variable Naming Conventions

### Standard Renames
- `fasta` → **`assembly`**: "Assembled contigs in FASTA format"
- `fastq`/`reads` → **`reads`**: "FASTQ reads (Illumina or Nanopore)"

### Channel Naming
- Use descriptive names: `ch_results`, `ch_logs`, `ch_nf_logs`, `ch_versions`
- For component-specific outputs: use the output type (e.g., `aln`, `csv`, `report`)

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
- [Troubleshooting](../reference/02-troubleshooting.md) - For common error solutions
