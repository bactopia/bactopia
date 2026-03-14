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
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, report, results, logs, nf_logs, versions)
 * - `report`: A tab-delimited report of hits, for full details please see [Abricate - Output](https://github.com/tseemann/abricate#output)
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
 * @citation prokka
 *
 * @note Uses EMPTY_* placeholder files for optional parameters
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input proteins
 * FASTA file of trusted proteins to first annotate from (Optional)
 *
 * @input prodigal_tf
 * Training file to use for gene prediction (Optional)
 *
 * @output record(meta, annotations, gff, gbk, fna, faa, ffn, sqn, fsa, tbl, txt, tsv, blastdb, results, logs, nf_logs, versions)
 * - `gff`: Annotation in GFF3 format, containing both sequences and annotations
 * - `gbk`: Annotation in GenBank format, containing both sequences and annotations
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
- **Path-workarounds**: Uses EMPTY_* files for optional inputs
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
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation mlst
 *
 * @modules csvtk_concat, mlst
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
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
include { gather              } from 'plugin/nf-bactopia'

workflow MLST {
    take:
    assembly: Channel<Record>
    db: Path

    main:
    MLST_MODULE(assembly, db)
    CSVTK_CONCAT(gather(MLST_MODULE.out, 'mlst', field: 'tsv'), 'tsv', 'tsv')

    emit:
    sample_outputs = MLST_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
```

**Key Characteristics**:
- **2-channel emit**: `sample_outputs` (module record passthrough) and `run_outputs` (aggregated)
- **gather() with field**: Extracts `tsv` field from module records for aggregation
- **Typed take block**: `assembly: Channel<Record>`, `db: Path`

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

### Input Patterns for Modules

The codebase uses Record-typed inputs with explicit named parameters:

**Assembly-based modules** (single file input):
```groovy
// Single assembly file
input:
(_meta: Map, assembly: Path): Record
```

**Read-based modules** (multi-read input with explicit slots):
```groovy
// Explicit positional slots for different read types
input:
(_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record
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
(_meta: Map, assembly: Path, meta_file: Path): Record
```

**Additional inputs** are declared on separate lines:
```groovy
input:
(_meta: Map, assembly: Path): Record
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

Subworkflows emit two channels -- module records pass through directly:

```groovy
emit:
sample_outputs = MODULE.out
run_outputs = CSVTK_CONCAT.out
```

### Workflow Branching Pattern

```groovy
ch_final_results = ch_results.branch{ meta, _file ->
    run: meta.scope == 'run'
    sample: meta.scope == 'sample'
}
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
                        _meta: [name: "GCF_000017085"],
                        assembly: file("${params.test_data_dir}/data/species/staphylococcus_aureus/genome/GCF_000017085.fna")
                    )
                )
                input[1] = file("${params.test_data_dir}/data/datasets/mlst/mlst.tar.gz")
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
