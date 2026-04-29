/**
 * Search, validate, gather, and standardize input samples.
 *
 * This subworkflow processes raw input samples through validation, standardization, and metadata
 * collection. It handles various input types including local FASTQ files, SRA/ENA accessions,
 * NCBI assembly accessions, and assemblies. The workflow can merge multiple sequencing runs,
 * download remote data, and simulate reads from assemblies using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art).
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1_files, r2_files, se_files, lr_files) with Set<Path> slots (pre-merge)
 * - Output: record(meta, r1, r2, se, lr) with Path? slots (post-merge, consolidated)
 *
 * @status stable
 * @keywords validation, download, merging, simulation, metadata, fastq, sra, ena, art
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation,resource-download,conditional-logic,no-test
 * @citation art, fastq_dl, fastq_scan, ncbigenomedownload, pigz
 *
 * @modules bactopia_gather, csvtk_concat
 *
 * @input record(meta, r1_files, r2_files, se_files, lr_files)
 * - `meta`: Groovy Record containing sample information
 * - `r1_files`: Illumina R1 read files (Set, elements may be null)
 * - `r2_files`: Illumina R2 read files (Set, elements may be null)
 * - `se_files`: Single-end read files (Set, elements may be null)
 * - `lr_files`: Long read files (ONT/PacBio) or assembly for simulation (Set, elements may be null)
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited metadata file describing the valid samples
 *
 * @output run_outputs
 * - `csv`: Aggregated metadata from all samples
 *
 * @output reads
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 * - `fna?`: Assembly file for assembly-based samples
 */
// bactopia-lint: ignore S013
nextflow.enable.types = true

include { GATHER as GATHER_MODULE } from '../../../modules/bactopia/gather/main'
include { CSVTK_CONCAT            } from '../../../modules/csvtk/concat/main'
include { filterWithData          } from 'plugin/nf-bactopia'
include { gatherCsvtk             } from 'plugin/nf-bactopia'

workflow GATHER {
    take:
    samples: Channel<Record>

    main:
    ch_gather = GATHER_MODULE(samples)
    ch_csvtk_concat = CSVTK_CONCAT(gatherCsvtk(ch_gather, 'tsv', [name: 'meta']), 'tsv', 'tsv')

    emit:
    // Downstream inputs
    reads = filterWithData(ch_gather, ['r1', 'r2', 'se', 'lr', 'fna'])

    // Published outputs
    sample_outputs = ch_gather
    run_outputs = ch_csvtk_concat
}
