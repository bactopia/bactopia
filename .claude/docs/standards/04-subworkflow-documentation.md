# Subworkflow Documentation Methodology

## Overview
This guide provides a comprehensive methodology for creating consistent and accurate GroovyDoc documentation for Bactopia subworkflows. It combines the style guidelines from [01-style-guide.md](01-style-guide.md) with the decision logic from [02-logic-rules.md](02-logic-rules.md) and implementation details from [03-technical-specs.md](03-technical-specs.md).

## 1. Methodology
To ensure accuracy and consistency, follow this step-by-step process for each subworkflow:

### 1.1 Analyze `main.nf` (Subworkflow)
1. **`take` block:** Identify input channels to determine `input-type` (1 Channel = `single`, 2+ Channels = `multiple`).
2. **`include` statements:** Identify dependencies to populate `@modules` or `@subworkflows`.
3. **`emit` block (Verification Step):**
   - Confirm that named outputs map correctly to underlying module outputs.
   - **Crucial:** Verify that the aggregated channels (`results`, `logs`, `nf_logs`, `versions`) are constructed to include the corresponding outputs from **every** module/subworkflow called in the `main` block.
4. **Tag Verification:** Ensure all applicable tags are included:
   - `aggregation` - if the workflow merges samples (e.g., csvtk_concat)
   - `database-dependent` - if any tool requires an external database
   - `conditional-logic` - if the workflow has conditional branches (if/else) based on boolean parameters
   - `archive-output` - if the workflow creates compressed archives
   - `resource-download` - if the workflow downloads external resources

### 1.2 Analyze `main.nf` (Included Modules)
- Read the module code to understand the content of the files listed in the subworkflow's `emit` block.

### 1.3 Consult `data/citations.yml`
- Retrieve the official Tool Name, Citation Key (`@citation`), and Source Code Repository URL.

## 2. Style & Formatting Standards

### 2.1 Header & Description
- **Format:** One-sentence summary followed by a detailed description.
- **Links:** Use Markdown `[Tool Name](https://repo)` format.
  - **Protocol:** Always use `https`.
  - **Target:** Link to the **Source Code Repository** (e.g., GitHub), not documentation aggregators.
- **Updates:** Note that **Panaroo** is now the default for pangenome analysis (not PIRATE).

### 2.2 Tags
- **`@status`:** Always `stable`.
- **`@keywords`:** Comma-separated list relevant to the tool/biology.
- **`@tags`:**
  - **`complexity`:** Determined by workflow characteristics:
    - `simple`: Linear Input → Output. The process takes an input, runs a single command/utility, and produces an output. Little to no branching logic (if/else). Usually relies on a single tool. Very fast, low resource footprint.
    - `moderate`: Analysis with Logic or Parameters. Includes conditional execution, multiple command steps within the script, or significant parameter handling. Often generates multiple distinct types of outputs. Wraps a standard bioinformatics tool that typically requires a specific database or reference file.
    - `complex`: Orchestration and Heavy Lifting. Involves calling other subworkflows, complex branching (e.g., switching tools based on flags), or massive parallelization/aggregation patterns. Often involves joining multiple distinct data streams. High resource usage or long runtimes; often acts as a "pipeline-within-a-pipeline".
  - **`input-type`:**
    - `single`: The `take` block defines exactly **1 Channel** (plus any number of parameters).
    - `multiple`: The `take` block defines **2 or more Channels**.
  - **`output-type`:** Usually `multiple` for subworkflows.
  - **`features`:** Include `aggregation` if the workflow merges samples (e.g., `csvtk_concat`).
- **Dependencies:**
  - Use **`@modules`** if `include` points to `../modules/...`.
  - Use **`@subworkflows`** if `include` points to `../subworkflows/...`.

### 2.3 Inputs

- **Format:** `@input tuple(meta, variable_name)`
- **Standard Renaming Rules:**
    - `fasta` → **`assembly`**: "Assembled contigs in FASTA format"
    - `fastq`/`reads` → **`reads`**: "FASTQ reads (Illumina or Nanopore)"
    - `meta` → "Groovy Map containing sample information"

#### Explicit Positional Read Slots

For subworkflows accepting reads, use explicit positional slots for clarity:

```groovy
@input tuple(meta, r1, r2, se, lr)
- `meta`: Groovy Map containing sample information
- `r1`: Illumina R1 reads (paired-end forward)
- `r2`: Illumina R2 reads (paired-end reverse)
- `se`: Single-end Illumina reads
- `lr`: Long reads (ONT/PacBio)
```

This pattern uses `Tuple<Map, Path?, Path?, Path?, Path?>` where each slot is optional, allowing the subworkflow to handle different read configurations.

**Subworkflows using this pattern**: ariba, kraken2, bracken, scrubber, teton

### 2.4 Outputs

- **Source of Truth:** Document all outputs listed in the subworkflow's `emit` block.
- **Vertical Alignment:** Pad variable names with spaces so all descriptions start at the same column index.
- **Mandatory Aggregates:** Always include these four at the end, explicitly stating they aggregate **all** underlying processes:
    - `results`    (Aggregated results channel containing all output files)
    - `logs`       (Aggregated logs channel containing all execution logs)
    - `nf_logs`    (Aggregated Nextflow execution scripts and logs for debugging from all processes)
    - `versions`   (Aggregated version information from all executed tools)

#### Special Output Channels

Some subworkflows emit additional "special" outputs designed for pass-through to downstream subworkflows:

```groovy
@output special_tsv    Intermediate output for downstream subworkflow consumption
```

**Example**: The scrubber subworkflow emits `special_tsv` which is consumed by the teton subworkflow.

### 2.5 The @note Tag

Use `@note` to document special requirements, caveats, or important information:

```groovy
@note Database can be automatically downloaded or provided as pre-existing tarball
```

Common uses:
- Database requirements
- Optional features
- Conditional behavior
- Important caveats

## 3. Information Sources
- **Primary:** Subworkflow `main.nf` (for logic/flow) and Module `main.nf` (for file details).
- **Secondary:** `data/citations.yml` (for official repository links).
- **Reference:** Bactopia documentation (internal) and Tool Repositories (external).

## 4. Examples of Successfully Updated Subworkflows

### 4.1 Pangenome (Moderate Complexity with Conditional Logic)
```groovy
/**
 * Perform pangenome analysis with optional core-genome phylogeny.
 *
 * This subworkflow creates a pangenome from GFF3 annotation files using one of three
 * tools: [Panaroo](https://github.com/gtonkinhill/panaroo) (default),
 * [PIRATE](https://github.com/SionBayliss/PIRATE), or
 * [Roary](https://github.com/sanger-pathogens/roary). It generates core-genome alignments
 * and gene presence/absence matrices, followed by SNP distance calculations using
 * [snp-dists](https://github.com/tseemann/snp-dists). The workflow conditionally executes
 * the selected pangenome tool based on Boolean parameters.
 *
 * @status stable
 * @keywords alignment, core-genome, pan-genome, phylogeny, comparative genomics
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, conditional-logic
 * @citation pirate, panaroo, roary, snpdists
 *
 * @subworkflows pirate, roary, panaroo, snpdists
 *
 * @input tuple(meta, gff)
 * - `meta`: Groovy Map containing sample information
 * - `gff`: Set of GFF3 annotation files from assembled genomes
 *
 * @input use_pirate
 * Boolean flag to use PIRATE for pangenome analysis
 *
 * @input use_roary
 * Boolean flag to use Roary for pangenome analysis
 *
 * @output aln          Core-genome alignment file containing genes present across all input genomes
 * @output csv          Gene presence/absence matrix showing which genes are present in each genome
 * @output results      Aggregate of all result files from pangenome analysis and SNP distances
 * @output logs         Aggregate of all log files from executed tools
 * @output nf_logs      Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Software version information from all executed tools
 */
```
**Key Features:**
- Uses `@subworkflows` since all includes point to other subworkflows
- Tags include `conditional-logic` for if/else branches and `aggregation` for merging results
- Documents Panaroo as the default tool (updated from PIRATE)

### 4.2 ARIBA (Moderate Complexity with Resource Download)
```groovy
/**
 * Rapidly identify genes by creating local assemblies from paired-end reads.
 *
 * This subworkflow uses [ARIBA](https://github.com/sanger-pathogens/ariba)
 * (Antimicrobial Resistance Identification By Assembly) to rapidly identify genes
 * in a database by creating local assemblies. It first downloads and prepares an ARIBA database,
 * then analyzes paired-end reads to identify genes, and finally aggregates results across all samples.
 *
 * @status stable
 * @keywords bacteria, reads, antimicrobial resistance, virulence, local assembly
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent, resource-download
 * @citation ariba
 *
 * @modules ariba_getref, ariba_run, csvtk_concat
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end reads in FASTQ format
 *
 * @input db
 * Database name for ARIBA analysis (e.g., ncbi, card, vfdb, resfinder, argannot)
 *
 * @output report         Per-sample tab-delimited reports of gene findings
 * @output summary        Per-sample CSV summaries of gene detection results
 * @output merged_report  Consolidated report containing gene findings from all samples
 * @output merged_summary Consolidated summary containing detection results from all samples
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions       Aggregated version information from all executed tools
 */
```
**Key Features:**
- Uses `@modules` since all includes point to modules
- Tags include `resource-download` and `database-dependent` for ARIBA database prep
- Shows vertical alignment in output descriptions

### 4.3 Teton (Complex Orchestration)
```groovy
/**
 * Perform taxonomic classification and estimate bacterial genome sizes.
 *
 * This subworkflow processes raw sequencing reads through a taxonomic classification
 * pipeline using [Kraken2](https://github.com/DerrickWood/kraken2) and [Bracken](https://github.com/jenniferlu717/Bracken)
 * to estimate bacterial genome sizes and separate bacterial from non-bacterial organisms.
 * It first removes host reads using the scrubber subworkflow, then classifies reads,
 * and finally creates sample sheets with genome size estimates for downstream Bactopia analysis.
 *
 * @status stable
 * @keywords metagenomics, taxonomy, classification, kraken, bracken, genome size
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, database-dependent, conditional-logic
 * @citation kraken2, bracken
 *
 * @subworkflows scrubber, bracken
 * @modules bactopia_samplesheet, csvtk_join, csvtk_concat
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: FASTQ reads (Illumina or Nanopore)
 *
 * @input db
 * Optional Kraken2 database path for taxonomic classification
 *
 * @input use_srascrubber
 * Boolean flag to use SRA scrubber for host read removal
 *
 * @output bacteria_tsv           Per-sample TSV files containing bacterial organisms and their properties
 * @output merged_bacteria_tsv    Consolidated TSV file of all bacterial organisms across samples
 * @output nonbacteria_tsv        Per-sample TSV files containing non-bacterial organisms
 * @output merged_nonbacteria_tsv Consolidated TSV file of all non-bacterial organisms across samples
 * @output sizemeup               Per-sample TSV files with genome size estimates
 * @output merged_sizemeup        Consolidated TSV file of genome size estimates across samples
 * @output report                 Joined TSV file combining scrubber and classification results
 * @output results                Aggregated results channel containing all output files
 * @output logs                   Aggregated logs channel containing all execution logs
 * @output nf_logs                Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions               Aggregated version information from all executed tools
 */
```
**Key Features:**
- Complexity: `complex` due to multiple subworkflow calls and data joining
- Mix of `@subworkflows` and `@modules` in dependencies
- Multiple Boolean parameters for conditional execution
- Shows proper handling of both optional database and conditional logic flags

### 4.4 SsuisSero (Simple with Aggregation)
```groovy
/**
 * Predict serotypes of Streptococcus suis from genome assemblies.
 *
 * This subworkflow uses [SsuisSero](https://github.com/jimmyliu1326/SsuisSero) to predict
 * serotypes of *Streptococcus suis* strains from genome assemblies based on the presence
 * of specific capsular genes. It processes each sample individually and aggregates the
 * results into a single consolidated report.
 *
 * @status stable
 * @keywords streptococcus suis, serotype, typing, prediction, capsular genes
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation ssuissero
 *
 * @modules ssuissero, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing serotype predictions
 * @output merged_tsv  Consolidated TSV file containing serotype predictions from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
```
**Key Features:**
- Input renamed from `fasta` to `assembly` per standard
- Complexity: `moderate` due to multiple processing steps
- Simple aggregation pattern with csvtk_concat

### 4.5 SCCmec (Multiple Output Types)
```groovy
/**
 * Identify SCCmec elements in Staphylococcus aureus genomes.
 *
 * This subworkflow uses [SCCmec](https://github.com/rpetit3/sccmec) to identify the
 * Staphylococcal Cassette Chromosome mec (SCCmec) element in *Staphylococcus aureus*
 * assemblies. It predicts the type based on the presence of specific *mec* and *ccr*
 * gene complexes, generating detailed BLAST results and typing information.
 *
 * @status stable
 * @keywords sccmec, staphylococcus aureus, mrsa, antimicrobial resistance, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation sccmec
 *
 * @modules sccmec, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv              Per-sample TSV files with SCCmec typing results
 * @output merged_tsv       Consolidated TSV file containing SCCmec typing from all samples
 * @output targets          Per-sample BLAST results for target sequences
 * @output target_details   Per-sample detailed results for target matches
 * @output regions          Per-sample BLAST results for SCCmec regions
 * @output regions_details  Per-sample detailed results for SCCmec region matches
 * @output results          Aggregated results channel containing all output files
 * @output logs             Aggregated logs channel containing all execution logs
 * @output nf_logs          Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions         Aggregated version information from all executed tools
 */
```
**Key Features:**
- Multiple output types from a single module
- Shows comprehensive documentation of all outputs
- Vertical alignment in output descriptions for readability

### 4.6 Assembler (Complex with Automatic Tool Selection)
```groovy
/**
 * Assemble bacterial genomes using automated assembler selection.
 *
 * This subworkflow automatically selects the optimal assembly strategy based on input read types:
 * - **Short Paired-End Reads:** Uses [Shovill](https://github.com/tseemann/shovill) (SKESA/SPAdes wrapper)
 * - **Short Single-End Reads:** Uses [Shovill](https://github.com/rpetit3/shovill) (SKESA/SPAdes wrapper)
 * - **Long Reads:** Uses [Dragonflye](https://github.com/rpetit3/dragonflye) (Flye/Miniasm wrapper)
 * - **Hybrid Assembly:** Uses [Unicycler](https://github.com/rrwick/Unicycler) or Dragonflye with short-read polishing
 *
 * The workflow performs individual assemblies per sample and aggregates assembly statistics
 * across all samples using [assembly-scan](https://github.com/rpetit3/assembly-scan) for comprehensive quality assessment.
 *
 * @status stable
 * @keywords bacteria, assembly, hybrid, shovill, dragonflye, unicycler, illumina, nanopore
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, conditional-logic, alternative-execution
 * @citation any2fasta, assembly-scan, bwa, dragonflye, flash, flye, medaka, megahit, miniasm, minimap2, nanoq, pigz, pilon, racon, rasusa, raven, samclip, samtools, shovill, shovill-se, skesa, spades, unicycler, velvet
 *
 * @modules bactopia_assembler, csvtk_concat
 *
 * @input tuple(meta, fq, extra)
 * - `meta`: Groovy Map containing sample information
 * - `fq`: Primary reads (Illumina paired-end or Nanopore)
 * - `extra`: Secondary reads for hybrid assembly or polishing (Optional)
 *
 * @output fna        Assembled contigs in FASTA format
 * @output fna_fq     Tuple containing assembly and primary reads for downstream analysis
 * @output tsv        Per-sample tab-delimited assembly statistics (N50, length, coverage)
 * @output merged_tsv Consolidated assembly statistics report across all samples
 * @output results    Aggregated results channel containing all output files
 * @output logs       Aggregated logs channel containing all execution logs
 * @output nf_logs    Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions   Aggregated version information from all executed tools
 */
```
**Key Features:**
- Automatic assembler selection based on input read type (not Boolean flags)
- Supports short reads (Shovill), long reads (Dragonflye), and hybrid assembly (Unicycler)
- Uses `alternative-execution` feature tag for multiple tool options
- Three-element input tuple with optional extra reads for hybrid assembly

## 5. What to Avoid
- **Do not** reference `module.config` or `schema.json`.
- **Do not** include parameter defaults or version numbers in the text.
- **Do not** count non-Channel parameters (Strings, Booleans, Paths) toward the `input-type` count.
- **Do not** use `http` links; ensure SSL/TLS is used (`https`).

## 6. Common Documentation Errors to Avoid

### 6.1 Copy-Paste Errors
When creating documentation for new subworkflows, ensure you update ALL fields. Common mistakes include:
- **Wrong short description**: Leaving the original tool's description (e.g., "Mass screening of contigs for antimicrobial resistance" when documenting a variant caller)
- **Wrong long description**: Referencing the wrong tool or its functionality
- **Wrong citation**: Using the template's citation key instead of the correct one from `data/citations.yml`
- **Generic output descriptions**: Single-word descriptions like "Vcf" or "Bam" instead of meaningful descriptions

**Bad Example** (copy-paste error):
```groovy
/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 * ...
 * @citation abricate
 * ...
 * @output vcf  Vcf
 */
```

**Correct Example** (properly documented snippy):

```groovy
/**
 * Call variants against a reference genome using Snippy.
 *
 * This subworkflow performs rapid haploid variant calling from bacterial sequence reads
 * using [Snippy](https://github.com/tseemann/snippy). It maps reads to a reference genome,
 * identifies SNPs and indels, and generates consensus sequences.
 * ...
 * @citation snippy
 * ...
 * @output vcf  Filtered variant calls in VCF format
 */
```

### 6.2 Legacy Comments
Remove any old-style comment headers before the GroovyDoc block:

```groovy
// BAD: Legacy comment that should be removed
//
// busco - Assembly completeness based on evolutionarily informed expectations
//

/**
 * Assess genome assembly completeness using BUSCO.
 ...
```

Should be:

```groovy
/**
 * Assess genome assembly completeness using BUSCO.
 ...
```

### 6.3 Citation Verification
**Always verify citations against `data/citations.yml`** before finalizing documentation:
1. Open `data/citations.yml`
2. Search for the tool name under the `tools:` section
3. Use the exact key as your `@citation` value
4. For multiple tools, use comma-separated format: `@citation tool1, tool2`

### 6.4 Formatting Issues
- **No extra blank lines** between the closing `*/` and `nextflow.preview.types = true`
- **Consistent indentation** (4 spaces) within the GroovyDoc block
- **Vertical alignment** of output descriptions at the same column

## See Also
- [Style Guide](01-style-guide.md) - For general GroovyDoc templates and formatting
- [Logic Rules](02-logic-rules.md) - For complexity classification and tag definitions
- [Technical Specifications](03-technical-specs.md) - For variable naming conventions and implementation details