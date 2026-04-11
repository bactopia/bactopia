# Module Documentation Methodology

## Overview
This guide provides a comprehensive methodology for creating consistent and accurate GroovyDoc documentation for Bactopia modules. It defines the standards for documenting individual tool implementations that form the foundation of the Bactopia pipeline.

## 1. Module Architecture

### 1.1 Module Structure
Bactopia modules are individual process definitions that execute specific bioinformatics tools. Each module:
- Wraps a single tool or closely related functionality
- Accepts standardized inputs using Record syntax (e.g., `(meta: Record, assembly: Path): Record`)
- Emits standardized outputs (including `logs`, `nf_logs`, `versions`)
- Uses static typing throughout
- Handles optional parameters via Path? workarounds

### 1.2 Standard Module Components
- **Header Block**: GroovyDoc with comprehensive metadata
- **Process Definition**: Named in UPPER_CASE (e.g., `PROKKA`, `MLST`)
- **Input Block**: Typed input channels and parameters
- **Output Block**: Typed output channels
- **Script Block**: Tool execution logic with proper variable handling

## 2. GroovyDoc Standard

### 2.1 Required Structure
```groovy
/**
 * <Single-sentence summary of the module's purpose>.
 *
 * <Detailed description explaining the tool, its methodology, and outputs.
 * Include the tool name in brackets [ToolName] with a link to the source repository.
 * Can span multiple lines and explain the biological or computational context.
 *
 * @status <stable|beta|deprecated>
 * @keywords <comma, separated, relevant keywords>
 * @tags complexity:<level> input-type:<type> output-type:<type> features:<features>
 * @citation <tool_name>
 *
 * @note <Optional: Special requirements, warnings, or notes>
 *
 * @input <input_channel>
 * <Input description>
 *
 * @input record(meta, <input_name>)
 * - `meta`: Groovy Record containing sample information
 * - `<input_name>`: <Description of the input files>
 *
 * @output record(meta, <field1>, <field2>, results, logs, nf_logs, versions)
 * - `<field1>`: Description of tool-specific output field
 * - `<field2>`: Description of tool-specific output field
 *
 * @results <Optional: supplemental or directory name, only if results has extra files>
 * - `<pattern>`: Description of published file not covered by named fields
 */
```

## 3. Tag Definitions and Usage

### 3.1 Complexity Levels

#### Simple
- **Definition**: Linear execution with minimal conditional logic
- **Characteristics**:
  - Single tool execution
  - Straightforward input → output mapping
  - Minimal parameter handling
- **Examples**: mlst, clermontyping, sistr

#### Moderate
- **Definition**: Multiple steps, basic conditional logic, or database dependency
- **Characteristics**:
  - Database download/extraction
  - File format detection and handling
  - Multiple output formats
- **Examples**: quast, fastani, emmtyper

#### Complex
- **Definition**: Multiple optional inputs, complex conditional logic, archive creation
- **Characteristics**:
  - Multiple optional Path? parameters
  - Complex workflow-dependent behavior
  - Archive or database creation
  - Advanced file handling
- **Examples**: prokka, bakta

### 3.2 Input-Type Classification

#### None
- **Definition**: No sample/data channels; only parameters
- **Pattern**: No record input block; may accept simple `Path` or `String` parameters
- **Use Case**: Utility modules for downloads, database setup, or internal maintenance
- **Examples**: wget, ariba/getref, bactopia/datasets, amrfinderplus/update

#### Single Input
- **Definition**: One primary data channel (plus parameters)
- **Pattern**: `(meta: Record, assembly: Path): Record` for assemblies
- **Pattern**: `(meta: Record, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record` for read-based modules
- **Examples**: Most modules that process a single data type

#### Multiple Inputs
- **Definition**: Multiple data channels or complex multi-element records
- **Pattern**: `(meta: Record, assembly: Path, meta_file: Path): Record` for multiple required files
- **Examples**: quast (assembly + meta_file)

### 3.3 Output-Type Classification

#### Single Output
- **Definition**: One primary output channel (plus standard outputs)
- **Use Case**: Tools that produce a single result file
- **Examples**: mlst (tsv)

#### Multiple Outputs
- **Definition**: Multiple distinct output channels
- **Use Case**: Tools producing multiple file formats or result types
- **Examples**: prokka (gff, gbk, fna, faa, etc.)

### 3.4 Feature Tags

**Important**: Feature values must be comma-separated WITHOUT spaces (e.g., `features:database-dependent,conditional-logic`, NOT `features:database-dependent, conditional-logic`)

#### Technical Features
- **database-dependent**: Requires external database
- **conditional-logic**: Contains if/else statements or complex conditional processing
- **archive-output**: Creates compressed archives (tar/zip)
- **compression**: Handles file compression/decompression
- **conditional-input**: Accepts optional `Path?` inputs
- **internet-access**: Requires internet connection for downloads
- **alternative-execution**: Multiple tool options (e.g., assembler can use shovill, dragonflye, unicycler)
- **resource-download**: Downloads external databases, datasets, or files

#### Processing Features
- **filtering**: Filters input data based on criteria
- **custom-outputs**: Non-standard output patterns

### 3.5 Common @note Patterns

The `@note` tag documents special requirements or context. Common patterns include:

#### Database Notes
- **Database Required**: External database must be provided by user
- **Database Bundled**: Tool includes its own database
- **Database Included**: Multiple databases bundled (e.g., abricate)

#### Internet/Resource Notes
- **Internet Required**: Process requires active internet connection

#### Utility Module Notes
- **Internal Utility**: Module used for internal pipeline operations
- **Internal Maintenance**: Module used for dataset building/updating

**Example usage**:
```groovy
* @note Database Bundled
* Kleborate bundles the required databases for species identification, MLST,
* and virulence/resistance gene detection.
```

## 4. Input Documentation Standards

### 4.1 Primary Data Inputs (Assembly)

```groovy
@input record(meta, assembly)
- `meta`: Groovy Record containing sample information
- `assembly`: Assembled contigs in FASTA format
```

### 4.2 Read Inputs (Explicit Positional Slots)

For modules accepting reads, use explicit positional slots with `?` suffix on optional fields:

```groovy
@input record(meta, r1?, r2?, se?, lr?)
- `meta`: Groovy Record containing sample information
- `r1?`: Illumina R1 reads (paired-end forward)
- `r2?`: Illumina R2 reads (paired-end reverse)
- `se?`: Single-end Illumina reads
- `lr?`: Long reads (ONT/PacBio)
```

The `?` suffix mirrors the `Path?` type in the take block, indicating these fields may be null.

### 4.3 Database Parameters
```groovy
@input db
Directory or compressed tarball containing the <tool> database
```

### 4.4 Optional Parameters

Use `?` suffix on the parameter name instead of "(Optional)" in the description:

```groovy
@input proteins?
FASTA file of trusted proteins to first annotate from

@input prodigal_tf?
Training file to use for gene prediction
```

## 5. Output Documentation Standards

### 5.1 Record Output Format

Modules use `@output record(...)` to document their outputs. The record line lists ALL fields from the actual `record()` output block, followed by description lines for tool-specific fields only.

```groovy
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited summary containing the Sample, Scheme, ST, and Allele IDs
```

### 5.2 Formatting Rules

1. **Single `@output record(...)` line** listing all fields from the process `record()` output block
2. **Tool-specific fields only** get ` * - ` description lines
3. **Skip standard fields** -- do NOT describe: `meta`, `results`, `logs`, `nf_logs`, `versions`
4. **Skip convenience bundles** -- e.g., `annotations` (fna+faa+gff bundle) is listed in the record but not described since individual fields already cover it
5. **Single line per field** -- each description must fit on one line, no wrapping
6. **Field names must match** the actual `record()` output block exactly
7. **Optional output fields** -- add `?` to fields that use `optional: true` in their `file()` call or conditional null passthrough
8. **Standard fields never get `?`** -- even if `logs` uses `optional: true` for file matching, it is not semantically optional

### 5.3 Standard Fields (Never Described)

These fields appear in most `record()` outputs but are not documented with description lines:

- `meta` -- Groovy Record containing sample information and output paths
- `results` -- List of result files for publishing
- `logs` -- Optional software execution logs
- `nf_logs` -- Nextflow execution scripts and logs
- `versions` -- YAML file with software versions

### 5.4 Output Patterns

#### Single Tool-Specific Output
```groovy
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited summary containing the results
```

#### Multiple Tool-Specific Outputs
```groovy
 * @output record(meta, gff, gbk, fna, faa, ffn, sqn, fsa, tbl, txt, tsv, blastdb, annotations, results, logs, nf_logs, versions)
 * - `gff`: Annotation in GFF3 format, containing both sequences and annotations
 * - `gbk`: Annotation in GenBank format, containing both sequences and annotations
 * - `fna`: Nucleotide FASTA file of the input contig sequences
 * - `faa`: Protein FASTA file of the translated CDS sequences
 * - `ffn`: Nucleotide FASTA file of all prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
 * - `sqn`: An ASN1 format "Sequin" file for submission to GenBank
 * - `fsa`: Nucleotide FASTA file of the input contig sequences (with adjusted headers)
 * - `tbl`: Feature table file for GenBank submission
 * - `txt`: Summary statistics of the annotation
 * - `tsv`: Tab-separated file of all annotated features
 * - `blastdb`: A compressed tar.gz archive of BLAST databases created from the input
```

#### No Tool-Specific Outputs
```groovy
 * @output record(meta, results, logs, nf_logs, versions)
```

#### Download/Utility Modules
```groovy
 * @output record(db, logs)
 * - `db`: A compressed tarball of the latest AMRFinder+ database
```

### 5.5 Published Results Documentation (@results)

Some modules publish additional files via the `results` list that are NOT named record fields (e.g., `supplemental/` directories, tool output directories). Use `@results` to document these so users know what files appear in their output directory.

**When to use:** Only when the `results` list contains `files()` patterns beyond what's already covered by named output fields. If every `files()` pattern in results corresponds to a named field, no `@results` is needed.

**Format:**
```groovy
 * @output record(meta, masked_aln, results, logs, nf_logs, versions)
 * - `masked_aln`: The input alignment with recombinant regions masked
 *
 * @results supplemental
 * - `*.recombination_predictions.gff.gz`: Predicted recombination regions in GFF format
 * - `*.per_branch_statistics.csv`: Statistics for each branch of the phylogeny
 * - `*.final_tree.tre`: Final recombination-free phylogenetic tree
```

**Rules:**
1. `@results` goes AFTER `@output` (with a blank ` *` line between them)
2. The header word describes the source: `supplemental`, `additional`, or a directory name (e.g., `panaroo/`, `pirate/`)
3. Each bullet describes a file pattern and what it contains
4. Use the same ` * - \`pattern\`: Description` format as @output field descriptions
5. Keep descriptions concise (one line per item)

**Common patterns:**
- `@results supplemental` -- for modules with a `supplemental/` directory
- `@results additional` -- for extra files not in a specific directory
- `@results <dirname>/` -- for modules publishing a tool's full output directory (e.g., `panaroo/`)

## 6. Implementation Patterns

### 6.1 Path? Parameter Handling

Optional `Path?` parameters are `null` when absent. Use null checks:

```groovy
def proteins_opt = proteins != null ? "--proteins ${proteins.getName()}" : ""
```

#### Conditional Database Handling
```groovy
def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
if [ "${is_tarball}" == "true" ]; then
    # Extract tarball
else
    # Use directory directly
fi
```

### 6.2 Meta Variable Construction

Bind `_meta = meta` at the top of the script block to preserve the upstream record, then rebuild `meta` as a new `Record` using the `record(...)` constructor. Records are immutable, so all fields must be set in a single expression (no post-construction mutation).

```groovy
script:
def _meta = meta
prefix = task.ext.prefix ?: "${_meta.name}"

// Create a new meta variable
meta = record(
    id: "${prefix}-${task.process}",
    name: prefix,
    scope: task.ext.scope,
    output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
    logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
    process_name: task.ext.process_name
)
```

### 6.3 Workflow-Dependent Behavior

Because `Record` fields cannot be mutated after construction, workflow-dependent values are expressed as inline ternaries inside the `record(...)` call rather than via conditional re-assignment.

```groovy
meta = record(
    id: "${prefix}-${task.process}",
    name: prefix,
    scope: task.ext.wf == "pangenome" ? "run" : task.ext.scope,
    output_dir: task.ext.wf == "pangenome" ? "prokka/${prefix}" : "${prefix}/main/annotator/prokka/",
    logs_dir: task.ext.wf == "pangenome" ? "prokka/${prefix}/logs" : "${prefix}/main/annotator/prokka/logs/",
    process_name: task.ext.process_name
)
```

### 6.4 Version Information

```bash
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tool_name: $( echo $(tool --version 2>&1) | sed 's/^.*pattern //' )
    database_version: $(cat database/VERSION 2>/dev/null || echo "unknown")
END_VERSIONS
```

### 6.5 Hardcoded Version Strings

Some tools do not provide version information via CLI. In these cases, use a hardcoded VERSION variable with a comment explaining why:

```groovy
// WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
def VERSION = '2.1'
```

**Modules using this pattern**: gamma, mcroni, plasmidfinder, pneumocat, midas, ssuissero

When updating container versions for these modules, remember to manually update the VERSION string.

### 6.6 Stage Blocks

Some modules use `stage:` blocks to organize input file staging. This is useful when inputs need specific directory structures:

```groovy
stage:
stageAs "input/*", assembly
```

**Multi-file staging example** (stecfinder):

```groovy
stage:
stageAs 'fna/*', fna
stageAs 'reads/r1/*', r1
stageAs 'reads/r2/*', r2
stageAs 'reads/se/*', se
stageAs 'reads/lr/*', lr
```

**Modules using stage blocks**: prokka, agrvate, fastani, roary, pirate, stecfinder

### 6.7 Runtime Meta Fields

Some modules add runtime-determined fields to the meta map for internal use:

```groovy
meta.is_paired = has_r1 && has_r2
meta.single_end = has_se && !has_r1 && !has_r2
meta.runtype = _meta.runtype
```

These fields are computed at runtime based on which inputs are provided.

## 7. Examples by Complexity

### 7.1 Simple Module (mlst)
```groovy
/**
 * Automatic Multi-Locus Sequence Typing (MLST) of genome assemblies.
 *
 * Uses [mlst](https://github.com/tseemann/mlst) to scan genome assemblies against traditional
 * PubMLST schemes. It automatically detects the likely species scheme, identifies the alleles
 * for the 7 housekeeping genes, and assigns a Sequence Type (ST).
 *
 * @status stable
 * @keywords bacteria, typing, mlst, sequence type, pubmlst, alleles
 * @tags complexity:simple input-type:single output-type:single features:database-dependent,conditional-logic
 * @citation mlst
 *
 * @note Database Required
 * Requires the MLST database (derived from PubMLST) to be available.
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * Directory or compressed tarball containing the MLST database schemes
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited summary containing the Sample, Scheme, ST, and Allele IDs
 */
```

### 7.2 Moderate Module (quast)
```groovy
/**
 * Quality Assessment Tool for Genome Assemblies.
 *
 * Uses [QUAST](https://github.com/ablab/quast) to evaluate genome assemblies by computing various
 * metrics such as N50, gene counts, and assembly length.
 *
 * @status stable
 * @keywords quast, assembly, quality control, n50, metrics
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation quast
 *
 * @input record(meta, assembly, meta_file)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 * - `meta_file`: Meta file containing reference size information
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Transposed report in TSV format
 *
 * @results supplemental
 * - `report.html`: Interactive HTML report with assembly metrics and plots
 * - `report.txt`: Plain text version of the assembly report
 * - `icarus.html`: Icarus contig size viewer
 * - `basic_stats/`: Directory containing GC content and coverage statistics
 */
```

### 7.3 Complex Module (prokka)
```groovy
/**
 * Annotate prokaryotic genomes.
 *
 * Uses [Prokka](https://github.com/tseemann/prokka) to rapidly annotate bacterial, archaeal,
 * and viral genomes, producing standards-compliant output files including GFF3, GenBank, and Sequin.
 *
 * @status stable
 * @keywords prokka, annotation, prokaryotic, bacteria, genbank, gff
 * @tags complexity:complex input-type:multiple output-type:multiple features:archive-output,compression,conditional-logic
 * @citation prokka
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input proteins?
 * FASTA file of trusted proteins to first annotate from
 *
 * @input prodigal_tf?
 * Training file to use for gene prediction
 *
 * @output record(meta, gff, gbk, fna, faa, ffn, sqn, fsa, tbl, txt, tsv, blastdb, annotations, results, logs, nf_logs, versions)
 * - `gff`: Annotation in GFF3 format, containing both sequences and annotations
 * - `gbk`: Annotation in GenBank format, containing both sequences and annotations
 * - `fna`: Nucleotide FASTA file of the input contig sequences
 * - `faa`: Protein FASTA file of the translated CDS sequences
 * - `ffn`: Nucleotide FASTA file of all prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
 * - `sqn`: An ASN1 format "Sequin" file for submission to GenBank
 * - `fsa`: Nucleotide FASTA file of the input contig sequences (with adjusted headers)
 * - `tbl`: Feature table file for GenBank submission
 * - `txt`: Summary statistics of the annotation
 * - `tsv`: Tab-separated file of all annotated features
 * - `blastdb`: A compressed tar.gz archive of BLAST databases created from the input
 */
```

## 8. Special Cases

### 8.1 Duplicate Output Fields
Some modules create duplicate record fields for pipeline routing:
```groovy
output:
record(
    meta: meta,
    classified: files("*.classified*.fastq.gz", optional: true),
    unclassified: files("*.unclassified*.fastq.gz", optional: true),
    classified_extra: files("*.classified*.fastq.gz", optional: true),
    unclassified_extra: files("*.unclassified*.fastq.gz", optional: true),
    results: [...],
    logs: files("*.{log,err}", optional: true),
    nf_logs: files(".command.*"),
    versions: files("versions.yml")
)
```

### 8.2 Workflow-Dependent Processing
```groovy
if (task.ext.wf == "scrubber" || task.ext.wf == "teton") {
    // Special processing for scrubbing workflow
}
```

### 8.3 Multi-Part Modules
Some tools are split across multiple modules (e.g., bakta/download, bakta/run):
- Each module should be fully documented
- Cross-reference related modules in @note if helpful

### 8.4 Utility/Setup Modules
Some modules are used for setup, downloading, or internal maintenance tasks rather than sample processing. These modules may have non-standard output structures:

**Examples**: `wget`, `ariba/getref`, `amrfinderplus/update`, `bakta/download`, `bactopia/datasets`

**Characteristics**:
- May not include `nf_logs` and `versions` as separate outputs (logs may be bundled in a subdirectory)
- Often have `input-type:none` or `input-type:single` (for parameter-only inputs)
- Typically include `features:resource-download,internet-access`
- Output structure focuses on the downloaded/created resource rather than sample results

**Example**:
```groovy
/**
 * Download and index the latest AMRFinder+ database.
 *
 * @status stable
 * @keywords bacteria, database, antimicrobial resistance, update, download, ncbi
 * @tags complexity:simple input-type:none output-type:single features:internet-access,archive-output,compression,database-dependent
 * @citation amrfinderplus
 *
 * @note Internal Maintenance
 * This process is primarily used internally by Bactopia to build and update the built-in datasets.
 *
 * @output record(db, logs)
 * - `db`: A compressed tarball of the latest AMRFinder+ database
 */
```

## 9. Schema.json Structure

Each module includes a `schema.json` file that defines parameter validation for the tool. This schema is used by Nextflow's schema validation system.

### 9.1 Template Structure

```json
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://raw.githubusercontent.com/bactopia/bactopia/master/modules/{tool_name}/schema.json",
    "title": "{Tool Name} Module",
    "description": "A module for {brief description of functionality}",
    "type": "object",
    "$defs": {
        "{tool}_parameters": {
            "title": "{Tool Name} Parameters",
            "type": "object",
            "description": "",
            "default": "",
            "fa_icon": "fas fa-exclamation-circle",
            "properties": {
                "{tool}_{param_name}": {
                    "type": "string|integer|boolean|number",
                    "default": "",
                    "description": "Description of the parameter",
                    "fa_icon": "fas fa-expand-arrows-alt"
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

### 9.2 Parameter Types

| JSON Type | Use For | Example |
|-----------|---------|---------|
| `"string"` | Text values, file paths | `"default": ""` or `"default": "ncbi"` |
| `"integer"` | Whole numbers | `"default": 80` |
| `"number"` | Decimal numbers | `"default": 0.95` |
| `"boolean"` | True/false flags | `"default": false` |

### 9.3 Naming Conventions

- **Parameter names**: Prefix with tool name (e.g., `mlst_minid`, `abricate_db`)
- **$defs key**: Use `{tool}_parameters` format
- **$id URL**: Point to raw GitHub path for the module's schema

### 9.4 Complete Example (mlst)

```json
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://raw.githubusercontent.com/bactopia/bactopia/master/modules/mlst/schema.json",
    "title": "MLST Module",
    "description": "A module for automatic MLST calling from assembled contigs",
    "type": "object",
    "$defs": {
        "mlst_parameters": {
            "title": "MLST Parameters",
            "type": "object",
            "description": "",
            "default": "",
            "fa_icon": "fas fa-exclamation-circle",
            "properties": {
                "mlst_scheme": {
                    "type": "string",
                    "default": "",
                    "description": "Don't autodetect, force this scheme on all inputs",
                    "fa_icon": "fas fa-expand-arrows-alt"
                },
                "mlst_minid": {
                    "type": "integer",
                    "default": 95,
                    "description": "Minimum DNA percent identity of full allele to consider 'similar'",
                    "fa_icon": "fas fa-expand-arrows-alt"
                },
                "mlst_mincov": {
                    "type": "integer",
                    "default": 10,
                    "description": "Minimum DNA percent coverage to report partial allele at all",
                    "fa_icon": "fas fa-expand-arrows-alt"
                },
                "mlst_minscore": {
                    "type": "integer",
                    "default": 50,
                    "description": "Minimum score out of 100 to match a scheme",
                    "fa_icon": "fas fa-expand-arrows-alt"
                },
                "mlst_nopath": {
                    "type": "boolean",
                    "default": false,
                    "description": "Strip filename paths from FILE column",
                    "fa_icon": "fas fa-expand-arrows-alt"
                },
                "mlst_db": {
                    "type": "string",
                    "default": "",
                    "description": "A custom MLST database to use, either a tarball or a directory",
                    "fa_icon": "fas fa-expand-arrows-alt"
                }
            }
        }
    },
    "allOf": [
        {
            "$ref": "#/$defs/mlst_parameters"
        }
    ]
}
```

### 9.5 Multi-Process Module Schemas

For modules with subdirectories (e.g., `abricate/run/`, `abricate/summary/`):

- Each subdirectory has its own `schema.json`
- Parameter names remain consistent with the tool prefix
- Shared parameters can be duplicated or split logically

### 9.6 Schema Location

```bash
modules/{tool_name}/schema.json           # Simple modules
modules/{tool_name}/{process}/schema.json # Multi-process modules
```

## 10. Module Configuration (module.config)

Each module includes a `module.config` file that defines parameter defaults, container images, and process-level settings via `task.ext` properties.

### 10.1 Template Structure

```groovy
params {
    // {tool}
    {tool}_flag = false
    {tool}_param1 = ""
    {tool}_param2 = 95
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
            params.{tool}_flag ? "--flag" : "",
            "--param2 ${params.{tool}_param2}"
        ].join(' ').replaceAll("\\s{2,}", " ").trim()

        // Environment information
        ext.toolName = "bioconda::{package}={version}".replace("=", "-").replace(":", "-").replace(" ", "-")
        ext.docker = "biocontainers/{package}:{version}--{build}"
        ext.image = "https://depot.galaxyproject.org/singularity/{package}:{version}--{build}"
        ext.condaDir = "${params.condadir}"
    }
}
```

**Section ordering within `withName` blocks:**

1. Routing — `ext.wf`, `ext.scope`, `ext.subdir`, `ext.logs_subdir`, `ext.process_name`
2. `// Tool arguments` — `ext.args` (and `ext.args2`, `ext.args3` if needed)
3. `// Environment information` — `ext.toolName`, `ext.docker`, `ext.image`, `ext.condaDir`
4. `// Module-specific parameters` — any additional `ext.*` properties (optional)

**Params block conventions:**

- Parameters must be in **alphabetical order**
- First line is a section comment using the module name in snake_case: `// {tool}` or `// {tool}_{process}`
- Modules with no parameters use `// No parameters` (capital N)
- `fa_icon` in schema.json is determined by type: `string` = `fas fa-font`, `integer` = `fas fa-hashtag`, `number` = `fas fa-percentage`, `boolean` = `fas fa-toggle-on`

### 10.2 Key Properties

| Property | Purpose | Example |
|----------|---------|---------|
| `ext.args` | CLI arguments passed to the tool | `"--minid 95 --mincov 10"` |
| `ext.toolName` | Conda package spec (used for conda env path) | `"bioconda-mlst-2.32.2"` |
| `ext.docker` | Docker container image reference | `"biocontainers/mlst:2.32.2--hdfd78af_0"` |
| `ext.image` | Singularity/Apptainer image URL | `"https://depot.galaxyproject.org/singularity/mlst:2.32.2--hdfd78af_0"` |
| `ext.condaDir` | Base directory for conda environments | `"${params.condadir}"` |
| `ext.wf` | Current workflow name | `params.wf` |
| `ext.scope` | Output scope: `"sample"` or `"run"` | `"sample"` |
| `ext.subdir` | Subdirectory under output path | `""` or `params.abricate_db` |
| `ext.logs_subdir` | Subdirectory for logs | `""` |
| `ext.process_name` | Tool name used in output directory paths | `"mlst"` |

### 10.3 No-Parameters Variant

When a module has no user-configurable parameters:

```groovy
params {
    // No parameters
}

process {
    withName: 'SSUISSERO' {
        ext.wf = params.wf
        ext.scope = "sample"
        ext.subdir = ""
        ext.logs_subdir = ""
        ext.process_name = "ssuissero"

        // Tool arguments
        ext.args = ""

        // Environment information
        ext.toolName = "bioconda::ssuissero=1.0.1".replace("=", "-").replace(":", "-").replace(" ", "-")
        ext.docker = "biocontainers/ssuissero:1.0.1--hdfd78af_1"
        ext.image = "https://depot.galaxyproject.org/singularity/ssuissero:1.0.1--hdfd78af_1"
        ext.condaDir = "${params.condadir}"

        // Version information not provided by tool on CLI
        ext.version = "1.0.1"
    }
}
```

### 10.4 Hardcoded Version Pattern

Some tools do not provide version information via CLI. In these cases, add `ext.version` to `module.config` and reference it in `main.nf`:

**module.config:**
```groovy
ext.version = "1.0.1"
```

**main.nf (script block):**
```groovy
// WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
def VERSION = '1.0.1'
```

### 10.5 Container URL Construction

Container URLs follow a deterministic pattern derived from the bioconda package:

- **Conda spec**: `bioconda::{package}={version}`
- **Docker**: `biocontainers/{package}:{version}--{build_hash}`
- **Singularity**: `https://depot.galaxyproject.org/singularity/{package}:{version}--{build_hash}`

The `{build_hash}` (e.g., `hdfd78af_0`) comes from the Anaconda API. Prefer `linux-64` builds; fall back to `noarch` if no `linux-64` build exists.

The `.replace("=", "-").replace(":", "-").replace(" ", "-")` chain on `ext.toolName` converts the conda spec into a valid directory name for the conda environment path.

## 11. Module Tests

Each module includes a test suite using the [nf-test](https://www.nf-test.com/) framework. For detailed testing guidance, see [Testing Framework](../project/04-testing-framework.md).

### 11.1 Test File Structure

```text
modules/{tool}/tests/
    main.nf.test       # Test specification
    main.nf.test.snap  # Snapshot file (generated by nf-test)
    nextflow.config    # Nextflow configuration for test
    nf-test.config     # nf-test runner configuration
```

For multi-process modules, each subdirectory has its own test suite:
```text
modules/{tool}/{process}/tests/
    main.nf.test
    nextflow.config    # Note: includeConfig paths are one level deeper
    nf-test.config
```

### 11.2 main.nf.test Template

```groovy
nextflow_process {
    name "Test {PROCESS_NAME}"
    script "../main.nf"
    process "{PROCESS_NAME}"
    tag "modules"
    tag "{tool}"

    test("{tool} - module - {sample_id}") {
        when {
            params {
                test_data_dir = System.getenv("BACTOPIA_TESTS") ?: ""
            }
            process {
                """
                input[0] = Channel.of(
                    record(
                        meta: [name: "{sample_id}"],
                        assembly: file("${params.test_data_dir}/data/species/{species}/genome/{sample_id}.fna.gz")
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
                    record.tsv,
                    record.versions
                ).match() }
            )
        }
    }
}
```

**Key patterns:**
- `input[0]` for the primary record input, `input[1]` for additional inputs (databases, etc.)
- Output accessed via `process.out[0][0]` to get the first record
- Snapshot only tool-specific fields + `meta` + `versions`
- **Never snapshot**: `results`, `logs`, `nf_logs` (unstable across runs)
- Default test species: Portiera (`portiera/genome/GCF_000292685.fna.gz`)

### 11.3 nextflow.config Template

```groovy
// Minimal config for module-level testing
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
includeConfig "../../../conf/base.config"
includeConfig "../../../conf/profiles.config"
```

**Path depth for multi-process modules:** Use `../../../../conf/` instead of `../../../conf/` since the module.config is one level deeper (e.g., `modules/bakta/run/tests/nextflow.config`).

### 11.4 nf-test.config Template

This file is identical across all modules:

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

### 11.5 Test Data Quick Reference

| Species | Path | When to use |
|---------|------|-------------|
| **Portiera** | `portiera/genome/GCF_000292685.fna.gz` | **Default for all assembly modules** |
| S. aureus | `staphylococcus_aureus/genome/GCF_000017085.fna` | Species-specific typing tools |
| N. gonorrhoeae | `neisseria_gonorrhoeae/genome/GCF_001047255.fna` | Has .gz variant |
| H. influenzae | `haemophilus_influenzae/genome/GCF_900478275.fna` | Has .gz variant |
| Portiera reads | `portiera/illumina/SRR2838702_R{1,2}.fastq.gz` | **Default for read-based modules** |

All test data is stored externally and referenced via the `BACTOPIA_TESTS` environment variable pointing to the [bactopia-tests](https://github.com/bactopia/bactopia-tests) repository.

## 12. Quality Checklist

Before completing module documentation, verify:

- [ ] All required tags are present (@status, @keywords, @tags, @citation)
- [ ] Input types are correctly documented as records with meta
- [ ] **Parameter names in code match documentation** (e.g., if `take:` has `assembly:`, doc should say `@input record(meta, assembly)`)
- [ ] `@output record(...)` lists all fields from the actual `record()` output block
- [ ] Tool-specific fields have ` * - ` description lines
- [ ] Standard fields (meta, results, logs, nf_logs, versions) are NOT described
- [ ] `@results` present if module publishes files beyond named record fields (e.g., supplemental/)
- [ ] Complexity level accurately reflects implementation
- [ ] Feature tags match actual implementation
- [ ] Tool name includes link to source repository
- [ ] Optional parameters are marked as "(Optional)"
- [ ] Database requirements are clearly stated
- [ ] Any Path? workarounds are noted in @note

## 13. What to Avoid

- **Do not** include parameter defaults or version numbers in the description
- **Do not** reference `module.config` or `schema.json` in documentation
- **Do not** use `http` links; ensure SSL/TLS is used (`https`)
- **Do not** document internal implementation details unless relevant to users
- **Do not** include command-line examples in the GroovyDoc

## See Also
- [Style Guide](01-style-guide.md) - For general GroovyDoc templates
- [Logic Rules](02-logic-rules.md) - For complexity classification logic
- [Technical Specifications](03-technical-specs.md) - For Path? workarounds and implementation details
- [Subworkflow Documentation](04-subworkflow-documentation.md) - For subworkflow-specific documentation standards
