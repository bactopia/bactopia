# Examples

## Overview
This section provides detailed examples of different component types in Bactopia, with annotations explaining key patterns and conventions.

## Module Examples

### Simple Module (Abricate)

**File**: `modules/abricate/run/main.nf`

```groovy
/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * Screens assemblies for antimicrobial resistance and virulence genes using
 * [Abricate](https://github.com/tseemann/abricate). It bundles several databases
 * including NCBI, CARD, ResFinder, PlasmidFinder, ARG-ANNOT, and VFDB.
 *
 * @status stable
 * @keywords bacteria, assembly, fasta, antimicrobial resistance, virulence, plasmid, mobile genetic elements
 * @tags complexity:simple input-type:single output-type:single features:database-dependent
 * @citation abricate
 *
 * @note Database Included
 * Abricate bundles multiple databases including NCBI, CARD, ResFinder, PlasmidFinder,
 * ARG-ANNOT, and VFDB.
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited report of hits, for full details please see [Abricate - Output](https://github.com/tseemann/abricate#output)
 */
```

**Key Characteristics**:
- **Simple complexity**: Single tool, linear execution
- **Single input**: One assembly file
- **Single output**: Report file plus standard generic fields
- **Database-dependent**: Bundles multiple databases

### Complex Module (Prokka)

**File**: `modules/prokka/main.nf`

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
 * @citation prokka, aragorn, barrnap, cdhit, hmmer, infernal, minced, nhmmer, prodigal, rnammer, signalp
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input proteins?
 * FASTA file of trusted proteins to first annotate from
 *
 * @input prodigal_tf?
 * Training file to use for gene prediction
 *
 * @output record(meta, gff, gbff, fna, faa, ffn, sqn, fsa, tbl, txt, tsv, blastdb, results, logs, nf_logs, versions)
 * - `gff`: Annotation in GFF3 format, containing both sequences and annotations
 * - `gbff`: Annotation in GenBank format, containing both sequences and annotations
 * - `fna`: Nucleotide FASTA file of the input contig sequences
 * - `faa`: Protein FASTA file of the translated CDS sequences
 * - `ffn`: Nucleotide FASTA file of all prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
 * - `sqn`: An ASN1 format "Sequin" file for submission to GenBank
 * - `fsa`: Nucleotide FASTA file of the input contig sequences, used by tbl2asn
 * - `tbl`: Feature Table file for NCBI submission
 * - `txt`: Summary statistics relating to the annotated features found
 * - `tsv`: Tab-separated file of all features (locus_tag, ftype, len_bp, gene, EC_number, COG, product)
 * - `blastdb`: A compressed tar.gz archive of BLAST+ databases of the contigs, genes, and proteins
 */
```

**Key Characteristics**:
- **Complex complexity**: Multiple conditional logic paths, many options
- **Multiple inputs**: Assembly plus optional proteins and training files
- **Multiple outputs**: Many output files in different formats
- **Conditional-input**: Uses `Path?` types for optional inputs
- **Compression**: Handles compressed outputs
- **Archive-output**: Creates compressed BLAST database archives

## Subworkflow Example

### MLST Subworkflow

**File**: `subworkflows/mlst/main.nf`

```groovy
/**
 * Determine multilocus sequence types (MLST) from bacterial assemblies.
 *
 * This subworkflow uses [mlst](https://github.com/tseemann/mlst) to scan assembled
 * contigs against PubMLST typing schemes and determine sequence types (STs). It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords mlst, sequence typing, pubmlst, bacteria
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation,database-dependent
 * @citation mlst
 *
 * @modules csvtk_concat, mlst
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Record containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * PubMLST database to use for MLST typing
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited summary containing the Sample, Scheme, ST, and Allele IDs
 *
 * @output run_outputs
 * - `csv`: A merged TSV file with mlst results from all samples
 */
nextflow.preview.types = true

include { MLST as MLST_MODULE } from '../../modules/mlst/main'
include { CSVTK_CONCAT        } from '../../modules/csvtk/concat/main'
include { gatherCsvtk         } from 'plugin/nf-bactopia'

workflow MLST {
    take:
    assembly: Channel<Record>
    db: Value<Path>

    main:
    ch_mlst = MLST_MODULE(assembly, db)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_mlst, 'tsv', [name: 'mlst']), 'tsv', 'tsv')

    emit:
    sample_outputs = ch_mlst
    run_outputs = ch_csvtk_concat
}
```

**Key Characteristics**:
- **2-channel emit**: `sample_outputs` (module record passthrough) and `run_outputs` (aggregated)
- **gatherCsvtk() for csvtk-friendly aggregation**: extracts the `tsv` field from each record into a single set keyed by `[name: 'mlst']`
- **Typed take block**: `assembly: Channel<Record>`, `db: Value<Path>`

## Entry Workflow Example

### Pangenome Bactopia Tool

**File**: `workflows/bactopia-tools/pangenome/main.nf`

```groovy
/**
 * Pangenome analysis with optional core-genome phylogeny.
 *
 * This Bactopia Tool creates a pangenome from GFF3 annotation files using one of three
 * tools: [Panaroo](https://github.com/gtonkinhill/panaroo) (default),
 * [PIRATE](https://github.com/SionBayliss/PIRATE), or
 * [Roary](https://github.com/sanger-pathogens/roary). It generates core-genome alignments
 * and gene presence/absence matrices, followed by SNP distance calculations.
 * You can supplement your pangenome with completed genomes using the --species or
 * --accessions parameters, which downloads genomes from RefSeq and annotates them with
 * Prokka. A phylogeny based on the core-genome alignment is created by IQ-Tree, with
 * optional recombination masking using ClonalFrameML. Finally, pan-genome wide
 * association studies can be conducted using Scoary.
 *
 * @status stable
 * @keywords alignment, core-genome, pan-genome, phylogeny, comparative genomics, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,aggregation,conditional-logic
 * @citation clonalframeml, iqtree, iqtree_modelfinder, iqtree_ufboot, ncbigenomedownload, panaroo, pirate, prokka, roary, scoary
 *
 * @subworkflows utils_bactopia-tools, pangenome, ncbigenomedownload, prokka, clonalframeml, iqtree, scoary
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input use_pirate
 * Use PIRATE as the pangenome tool instead of Panaroo
 *
 * @input use_roary
 * Use Roary as the pangenome tool instead of Panaroo
 *
 * @input species
 * Species name used to supplement the pangenome with RefSeq assemblies
 *
 * @input accession
 * Single NCBI Assembly RefSeq accession to supplement the pangenome
 *
 * @input accessions
 * Path to a file listing NCBI Assembly accessions to supplement the pangenome
 *
 * @input prokka_proteins
 * Path to trusted protein FASTA for Prokka homology-based annotation of downloaded genomes
 *
 * @input prokka_prodigal_tf
 * Path to a Prodigal training file for Prokka gene prediction
 *
 * @input skip_phylogeny
 * Skip core-genome phylogeny construction with IQ-Tree
 *
 * @input skip_recombination
 * Skip recombination masking with ClonalFrameML
 *
 * @input scoary_traits
 * Path to a Scoary trait file to run pan-GWAS against the pangenome
 *
 * @section Pangenome Results
 * @publish *.aln                 Core-genome alignment file containing genes present across all input genomes
 * @publish *.csv                 Gene presence/absence matrix showing which genes are present in each genome
 * @publish *.tsv                 SNP distance matrix between all samples
 *
 * @section Phylogeny Results
 * @note Only created if --skip_phylogeny is not enabled
 * @publish *.treefile            Maximum likelihood phylogenetic tree in Newick format
 * @publish *.iqtree              IQ-Tree analysis report with model selection and support values
 * @publish *.log                 IQ-Tree execution log
 *
 * @section Recombination Analysis
 * @note Only created if --skip_recombination is not enabled
 * @publish *.masked.aln          Core-genome alignment with recombination regions masked
 *
 * @section Association Analysis
 * @note Only created if --scoary_traits is specified
 * @publish scoary/*              Scoary association analysis results and plots
 *
 * @section Panaroo Results
 * @note Only created when Panaroo is selected as the pangenome tool
 * @publish panaroo/*             Panaroo-specific output files including graph and statistics
 *
 * @section PIRATE Results
 * @note Only created when PIRATE is selected as the pangenome tool
 * @publish pirate/*              PIRATE-specific output files including gene families and clusters
 *
 * @section Roary Results
 * @note Only created when Roary is selected as the pangenome tool
 * @publish roary/*               Roary-specific output files including gene presence/absence matrices
 *
 * @section Execution Logs
 * @publish logs/pangenome/*      Pangenome tool execution logs (stdout/stderr)
 * @publish logs/clonalframeml/*  ClonalFrameML execution logs (if executed)
 * @publish logs/iqtree/*         IQ-Tree execution logs (if executed)
 * @publish logs/scoary/*         Scoary execution logs (if executed)
 * @publish logs/nf-*             Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml          Software version information
 */
```

**Key Characteristics**:
- **User-facing**: Published results for end users
- **Multiple sections**: Pangenome / Phylogeny / Recombination / Association / per-tool (Panaroo, PIRATE, Roary) / Execution Logs / Versions
- **Conditional outputs**: `@note` lines flag sections that only materialize under specific parameters (`--skip_phylogeny` off, `--skip_recombination` off, `--scoary_traits` set, tool selection)
- **Uses @publish**: Not @output for user-facing files

## Common Patterns

### Input Patterns for Modules

The codebase uses Record-typed inputs with explicit named parameters:

**Assembly-based modules** (single file input):
```groovy
// Single assembly file
input:
record (
    meta: Record,
    fna: Path
)
```

**Read-based modules** (multi-read input with explicit slots):
```groovy
// Explicit positional slots for different read types
input:
record (
    meta: Record,
    r1: Path?,
    r2: Path?,
    se: Path?,
    lr: Path?
)
```

Where each parameter represents:
- `r1`: Illumina R1 reads (paired-end forward)
- `r2`: Illumina R2 reads (paired-end reverse)
- `se`: Single-end Illumina reads
- `lr`: Long reads (ONT/PacBio)

**Multiple distinct inputs** (e.g., assembly + metadata):
```groovy
// Two required files
input:
record (
    meta: Record,
    fna: Path,
    meta_file: Path
)
```

**Additional inputs** are declared on separate lines:
```groovy
input:
record (
    meta: Record,
    fna: Path
)
db         : Path
proteins   : Path?
prodigal_tf: Path?
```

### Standard Output Pattern for Modules

Modules emit a single `record()` with named fields (for downstream access) and generic fields (for publishing):

```groovy
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

**Note**: Use `file()` for named fields (returns `Path`), `files()` for generic fields and wildcards (returns `Set<Path>`).

### Subworkflow Emit Pattern

Subworkflows emit two channels. Capture each process's output channel at the call site -- direct `.out` access is not used -- and pass the handles through the `emit:` block:

```groovy
main:
ch_module = MODULE(input)
ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_module, 'tsv', [name: 'module']), 'tsv', 'tsv')

emit:
sample_outputs = ch_module
run_outputs = ch_csvtk_concat
```

## Test Patterns

### Module Test Example

```groovy
nextflow_process {
    name "Test MLST"
    script "../main.nf"
    process "MLST"
    tag "modules"
    tag "mlst"

    test("mlst - module - GCF_000017085") {
        when {
            params {
                test_data_dir = System.getenv("BACTOPIA_TESTS") ?: ""
            }
            process {
                """
                input[0] = Channel.of(
                    record(
                        meta: [name: "GCF_000017085"],
                        fna: file("${params.test_data_dir}/species/staphylococcus_aureus/uncompressed/GCF_000017085/main/assembler/GCF_000017085.fna")
                    )
                )
                input[1] = file("${params.test_data_dir}/datasets/mlst/mlst.tar.gz")
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

**Key Characteristics**:
- Uses `nextflow_process` block for module-level tests
- Input uses `record()` syntax matching the module's input declaration
- Output accessed via `process.out[0][0]` to get the record
- Assertions use `record.fieldName` to access specific fields
- Snapshot matching with MD5 checksums for reproducibility

## Meta Record Examples

These show the `scope`, `output_dir`, and `logs_dir` field values that get passed into `meta = record(...)` for each scope.

### Sample Scope
```groovy
meta = record(
    // ...
    scope: 'sample',
    output_dir: "${prefix}/tools/tool-name/subdir",
    logs_dir: "${prefix}/tools/tool-name/subdir/logs/logs_subdir",
    // ...
)
```

### Run Scope
```groovy
meta = record(
    // ...
    scope: 'run',
    output_dir: "tool-name",
    logs_dir: "tool-name/logs/logs_subdir/subdir",
    // ...
)
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
