# Examples

## Overview
This section provides detailed examples of different component types in Bactopia, with annotations explaining key patterns and conventions.

## Module Examples

### Simple Module (Abricate)

**File**: `modules/abricate/run/main.nf`

```groovy
/**
 * Screen assemblies for antimicrobial resistance genes.
 *
 * This process screens assembled contigs against resistance gene databases
 * using [Abricate](https://github.com/tseemann/abricate). It quickly identifies
 * antimicrobial resistance, virulence genes, and other relevant sequences in
 * bacterial genomes.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, screening
 * @tags complexity:simple input-type:single output-type:multiple features:database-dependent
 * @citation abricate
 *
 * @input tuple(meta, assembly)
 * Input assembly and metadata. The meta map contains sample information
 * such as sample ID, name, and other metadata. The assembly set contains
 * the assembled contigs in FASTA format to be screened.
 *
 * @output report    Abricate screening results (TSV format with resistance gene hits)
 * @output logs      Optional tool log files containing execution information
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  Software version information in YAML format
 */
```

**Key Characteristics**:
- **Simple complexity**: Single tool, linear execution
- **Single input**: One assembly file
- **Multiple outputs**: Report file plus standard channels
- **Database-dependent**: Requires database selection

### Complex Module (Prokka)

**File**: `modules/prokka/run/main.nf`

```groovy
/**
 * Annotate bacterial genomes with functional information.
 *
 * This process annotates bacterial contigs or complete genomes using [Prokka](https://github.com/tseemann/prokka).
 * It rapidly calls genes, translates them, and searches them against multiple protein databases
 * to produce comprehensive annotation in various standard formats.
 *
 * @status stable
 * @keywords bacteria, annotation, genome, prokaryote, functional annotation
 * @tags complexity:complex input-type:single output-type:multiple features:path-workarounds conditional-input conditional-logic compression archive-output
 * @citation prokka
 *
 * @note Optional: Proteins and training files improve annotation quality
 *
 * @input tuple(meta, assembly)
 * Input assembly for annotation. The meta map contains sample information,
 * and the assembly set contains the assembled contigs in FASTA format.
 *
 * @input proteins
 * Optional protein sequences for homology search. When provided, these trusted
 * protein sequences are used to improve annotation accuracy through homology.
 *
 * @input prodigal_tf
 * Optional Prodigal training file. Species-specific training data that improves
 * gene prediction accuracy.
 *
 * @output annotations Complete annotation package containing GFF, GBK, FASTA, and other formats
 * @output gff         Genome annotation in GFF3 format (standard for genome browsers)
 * @output gbk         Genome annotation in GenBank format (rich format suitable for NCBI submission)
 * @output fna         Nucleotide sequences of all predicted genes
 * @output faa         Protein sequences of all predicted genes
 * @output gff         Protein sequences of all predicted genes
 * @output ffn         Nucleotide sequences of all predicted genes (with frames)
 * @output sqn         Sequin file for NCBI submission
 * @output fsa         Translated protein sequences for BLAST
 * @output tbl         Feature table file
 * @output txt         Prokka annotation summary
 * @output tsv         Tabular summary of annotations
 * @output err         Error log from Prokka
 * @output log         Execution log from Prokka
 * @output logs        Optional tool log files containing execution information
 * @output nf_logs     Nextflow execution scripts and logs for debugging
 * @output versions    Software version information in YAML format
 */
```

**Key Characteristics**:
- **Complex complexity**: Multiple conditional logic paths, many options
- **Single input**: One assembly file
- **Multiple outputs**: Many output files in different formats
- **Path-workarounds**: Uses EMPTY_* files for optional inputs
- **Conditional-input**: Accepts optional protein and training files
- **Compression**: Handles compressed outputs
- **Archive-output**: Creates compressed archives

## Subworkflow Example

### Pangenome Subworkflow

**File**: `subworkflows/pangenome/main.nf`

```groovy
/**
 * Perform pangenome analysis with optional core-genome phylogeny.
 *
 * This subworkflow creates a pangenome from GFF3 annotation files using one of three
 * tools: PIRATE (default), Roary, or Panaroo. It generates core-genome alignments
 * and gene presence/absence matrices, followed by SNP distance calculations using
 * snp-dists. The workflow conditionally executes the selected pangenome tool based
 * on Boolean parameters.
 *
 * @status stable
 * @keywords alignment, core-genome, pan-genome, phylogeny, comparative genomics
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation conditional-logic
 * @citation pirate, panaroo, roary, snpdists
 *
 * @subworkflows pirate, panaroo, roary, snpdists
 *
 * @note Optional: SNP distance calculation is always performed on core-genome alignment
 *
 * @input gff
 * GFF3 annotation files from assembled genomes. Each tuple contains metadata
 * about the sample and a set of GFF3 files representing the annotation.
 *
 * @input use_pirate
 * Use PIRATE for pangenome analysis (default). When true, executes PIRATE
 * which is suitable for highly diverse datasets.
 *
 * @input use_roary
 * Use Roary for pangenome analysis. When true, executes Roary which is
 * suitable for closely related bacterial genomes.
 *
 * @output aln          Core-genome alignment file containing genes present across all input genomes
 * @output csv          Gene presence/absence matrix showing which genes are present in each genome
 * @output results      Aggregate of all result files from pangenome analysis and SNP distances
 * @output logs         Aggregate of all log files from executed tools
 * @output nf_logs      Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Software version information from all executed tools
 */
```

**Key Characteristics**:
- **Moderate complexity**: Conditional logic between tools
- **Single input**: GFF annotation files
- **Multiple outputs**: Core alignment, matrix, and standard channels
- **Aggregation**: Combines results from multiple tools
- **Conditional-logic**: Chooses tool based on parameters

## Entry Workflow Example

### Pangenome Bactopia Tool

**File**: `workflows/bactopia-tools/pangenome/main.nf`

```groovy
/**
 * Bactopia Tool: Pangenome Analysis.
 *
 * Performs comprehensive pangenome analysis from bacterial genomes.
 * Creates gene presence/absence matrices and builds phylogenetic trees.
 *
 * @status stable
 * @keywords pangenome, comparative genomics, phylogeny
 * @citation pirate, panaroo, roary, iqtree, scoary, clonalframeml
 *
 * @subworkflows ncbigenomedownload, prokka, pangenome, clonalframeml, iqtree, scoary
 *
 * @note Optional: Requires trait file for GWAS analysis with SCOARY
 *
 * @input scoary_traits
 * Path to trait file for genome-wide association studies.
 *
 * @section Pangenome Analysis
 * @publish *.csv         Gene presence/absence matrix
 * @publish *.aln         Core-genome alignment
 * @publish *.tree        Phylogenetic tree file
 *
 * @section GWAS Analysis
 * @note Only created if --scoary_traits is specified
 * @publish scoary/*.csv  Association results between genes and traits
 * @publish scoary/*.txt  Statistical summary of GWAS results
 *
 * @section Recombination Analysis
 * @note Only created if --skip_recombination is false
 * @publish *.masked.aln Recombination-masked alignment
 * @publish clonalframe/*.json ClonalFrameML analysis results
 *
 * @section Versions
 * @publish versions.yml   Software version information
 */
```

**Key Characteristics**:
- **User-facing**: Published results for end users
- **Multiple sections**: Groups related outputs
- **Conditional outputs**: Some files only created with specific parameters
- **Uses @publish**: Not @output for user-facing files

## Common Patterns

### Input Pattern for Modules
```groovy
input:
(_meta, assembly): Tuple<Map, Set<Path>>
```

### Standard Output Pattern for Modules
```groovy
output:
logs        = tuple(meta, files("*.{log,err}", optional: true))
nf_logs     = tuple(meta, files(".command.*"))
versions    = tuple(meta, file("versions.yml"))
```

### Subworkflow Aggregation Pattern
```groovy
emit:
results    = flattenPaths([ch_results])
logs       = flattenPaths([ch_logs])
nf_logs    = flattenPaths([ch_nf_logs])
versions   = flattenPaths([ch_versions])
```

### Workflow Branching Pattern
```groovy
ch_final_results = ch_results.branch{ meta, _file ->
    run: meta.scope == 'run'
    sample: meta.scope == 'sample'
}
```

## Meta Map Examples

### Sample Scope
```groovy
meta.scope = 'sample'
meta.output_dir = "${prefix}/tools/tool-name/subdir"
meta.logs_dir = "${prefix}/tools/tool-name/subdir/logs/logs_subdir"
```

### Run Scope
```groovy
meta.scope = 'run'
meta.output_dir = "tool-name"
meta.logs_dir = "tool-name/logs/logs_subdir/subdir"
```

## Version Tracking Example
```groovy
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
END_VERSIONS
```

Expected output:
```yaml
"ABRICATE_RUN":
    abricate: 1.0.1
```

## See Also
- [Style Guide](../standards/01-style-guide.md) - For template formats
- [Logic Rules](../standards/02-logic-rules.md) - For classification criteria
- [Technical Specifications](../standards/03-technical-specs.md) - For implementation details