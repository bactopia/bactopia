# Module Documentation Methodology

## Overview
This guide provides a comprehensive methodology for creating consistent and accurate GroovyDoc documentation for Bactopia modules. It defines the standards for documenting individual tool implementations that form the foundation of the Bactopia pipeline.

## 1. Module Architecture

### 1.1 Module Structure
Bactopia modules are individual process definitions that execute specific bioinformatics tools. Each module:
- Wraps a single tool or closely related functionality
- Accepts standardized inputs (usually `Tuple<Map, Set<Path>>`)
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
 * @input tuple(meta, <input_name>)
 * - `meta`: Groovy Map containing sample information
 * - `<input_name>`: <Description of the input files>
 *
 * @output record(meta, <field1>, <field2>, results, logs, nf_logs, versions)
 * - `<field1>`: Description of tool-specific output field
 * - `<field2>`: Description of tool-specific output field
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
- **Pattern**: `String`, `Map`, or `Path` parameters without `Tuple<Map, ...>` channels
- **Use Case**: Utility modules for downloads, database setup, or internal maintenance
- **Examples**: wget, ariba/getref, bactopia/datasets, amrfinderplus/update

#### Single Input
- **Definition**: One primary data channel (plus parameters)
- **Pattern**: `Tuple<Map, Path>` for assemblies
- **Pattern**: `Tuple<Map, Path?, Path?, Path?, Path?>` for read-based modules (r1, r2, se, lr)
- **Examples**: Most modules that process a single data type

#### Multiple Inputs
- **Definition**: Multiple data channels or complex multi-element tuples
- **Pattern**: `Tuple<Map, Path, Path>` for multiple required files
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
- **path-workarounds**: Uses EMPTY_* files for optional parameters
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

#### Technical Notes
- **Uses EMPTY_* placeholder files for optional parameters**: Documents Path? workarounds

**Example usage**:
```groovy
* @note Database Bundled
* Kleborate bundles the required databases for species identification, MLST,
* and virulence/resistance gene detection.
```

## 4. Input Documentation Standards

### 4.1 Primary Data Inputs (Assembly)

```groovy
@input tuple(meta, assembly)
- `meta`: Groovy Map containing sample information
- `assembly`: Assembled contigs in FASTA format
```

### 4.2 Read Inputs (Explicit Positional Slots)

For modules accepting reads, use explicit positional slots:

```groovy
@input tuple(meta, r1, r2, se, lr)
- `meta`: Groovy Map containing sample information
- `r1`: Illumina R1 reads (paired-end forward)
- `r2`: Illumina R2 reads (paired-end reverse)
- `se`: Single-end Illumina reads
- `lr`: Long reads (ONT/PacBio)
```

This pattern provides clear documentation of which read types are supported and uses `Path?` types for optional slots.

### 4.3 Database Parameters
```groovy
@input db
Directory or compressed tarball containing the <tool> database
```

### 4.4 Optional Parameters
```groovy
@input proteins
FASTA file of trusted proteins to first annotate from (Optional)

@input prodigal_tf
Training file to use for gene prediction (Optional)
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

### 5.3 Standard Fields (Never Described)

These fields appear in most `record()` outputs but are not documented with description lines:

- `meta` -- Groovy Map containing sample information and output paths
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

## 6. Implementation Patterns

### 6.1 Path? Parameter Handling

Two main approaches for optional parameters:

#### EMPTY_* File Detection (Preferred)

```groovy
def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ?
    "--proteins ${proteins.getName()}" : ""
```

**Note**: Use `.getName()` directly on `Path?` parameters. The older `.toList()[0].getName()` pattern is deprecated.

#### Conditional Database Handling
```groovy
def is_tarball = db.toList()[0].getName().endsWith(".tar.gz") ? true : false
if [ "${is_tarball}" == "true" ]; then
    # Extract tarball
else
    # Use directory directly
fi
```

### 6.2 Meta Variable Construction
```groovy
meta = [:]
meta.id = "${prefix}-${task.process}"
meta.name = prefix
meta.scope = task.ext.scope
meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
meta.process_name = task.ext.process_name
```

### 6.3 Workflow-Dependent Behavior
```groovy
if (task.ext.wf == "pangenome") {
    meta.scope = "run"
    meta.output_dir = "prokka/${prefix}"
}
else {
    meta.scope = "sample"
    meta.output_dir = "${prefix}/main/annotator/prokka/"
}
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
def VERSION = '2.1'
// WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
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
 * @input tuple(meta, assembly, meta_file)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 * - `meta_file`: Meta file containing reference size information
 *
 * @output record(meta, tsv, supplemental, results, logs, nf_logs, versions)
 * - `tsv`: Transposed report in TSV format
 * - `supplemental`: Supplemental files including plots and HTML reports
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
 * @note Uses EMPTY_* placeholder files for optional parameters
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input proteins
 * FASTA file of trusted proteins to first annotate from (Optional)
 *
 * @input prodigal_tf
 * Training file to use for gene prediction (Optional)
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

### 8.1 Duplicate Output Channels
Some modules create duplicate channels for pipeline routing:
```groovy
output:
classified     = tuple(meta, files("*.fastq.gz"))
classified_extra = tuple(meta, files("*.fastq.gz"), file("EMPTY_EXTRA", optional: true))
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

## 10. Quality Checklist

Before completing module documentation, verify:

- [ ] All required tags are present (@status, @keywords, @tags, @citation)
- [ ] Input types are correctly documented as tuples with meta
- [ ] **Parameter names in code match documentation** (e.g., if `take:` has `assembly:`, doc should say `@input tuple(meta, assembly)`)
- [ ] `@output record(...)` lists all fields from the actual `record()` output block
- [ ] Tool-specific fields have ` * - ` description lines
- [ ] Standard fields (meta, results, logs, nf_logs, versions) are NOT described
- [ ] Complexity level accurately reflects implementation
- [ ] Feature tags match actual implementation
- [ ] Tool name includes link to source repository
- [ ] Optional parameters are marked as "(Optional)"
- [ ] Database requirements are clearly stated
- [ ] Any Path? workarounds are noted in @note

## 11. What to Avoid

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
