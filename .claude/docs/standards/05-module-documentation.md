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
 * @output <output_name>    <Description of the output>
 * @output logs             Optional software execution logs containing warnings/errors
 * @output nf_logs          Nextflow execution scripts and logs for debugging
 * @output versions         A YAML formatted file with software versions
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

### 5.1 Tool-Specific Outputs
- **Name**: Clear, descriptive channel name
- **Description**: Concise explanation of what the output contains
- **Format**: Mention file format if not obvious from name

### 5.2 Standard Outputs (All Modules)
```groovy
@output logs         Optional software execution logs containing warnings/errors
@output nf_logs      Nextflow execution scripts and logs for debugging
@output versions     A YAML formatted file with software versions
```

### 5.3 Output Patterns

#### Single File Outputs
```groovy
@output tsv      A tab-delimited summary containing the results
```

#### Multiple File Collections
```groovy
@output supplemental  Supplemental files including plots and HTML reports
```

#### Archives
```groovy
@output blastdb    A compressed tar.gz archive of BLAST databases
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
 * @output tsv       A tab-delimited summary containing the Sample, Scheme, ST, and Allele IDs
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
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
 * @output tsv          Transposed report in TSV format
 * @output supplemental Supplemental files including plots and HTML reports
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
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
 * @output annotations  A tuple containing the FASTA, protein FASTA, and GFF3 files
 * @output gff          Annotation in GFF3 format, containing both sequences and annotations
 * @output gbk          Annotation in GenBank format, containing both sequences and annotations
 * @output fna          Nucleotide FASTA file of the input contig sequences
 * @output faa          Protein FASTA file of the annotated genes
 * @output ffn          Nucleotide FASTA file of the annotated genes
 * @output sqn          An ASN1 format "Sequin" file for submission to GenBank
 * @output fsa          Nucleotide FASTA file of the input contig sequences (with adjusted headers)
 * @output tbl          Feature table file
 * @output txt          Summary statistics about the annotation
 * @output tsv          Tab-separated file of all features
 * @output blastdb      A compressed tar.gz archive of BLAST databases created from the input
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
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
 * @output db       A compressed tarball of the latest AMRFinder+ database
 * @output logs     Optional software execution logs containing warnings/errors
 */
```

## 9. Quality Checklist

Before completing module documentation, verify:

- [ ] All required tags are present (@status, @keywords, @tags, @citation)
- [ ] Input types are correctly documented as tuples with meta
- [ ] All outputs are documented with clear descriptions
- [ ] Standard outputs (logs, nf_logs, versions) are included
- [ ] Complexity level accurately reflects implementation
- [ ] Feature tags match actual implementation
- [ ] Tool name includes link to source repository
- [ ] Optional parameters are marked as "(Optional)"
- [ ] Database requirements are clearly stated
- [ ] Any Path? workarounds are noted in @note

## 10. What to Avoid

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
