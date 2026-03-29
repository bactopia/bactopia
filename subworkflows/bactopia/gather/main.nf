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
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, resource-download, conditional-logic, no-test
 * @citation art, fastq-dl, fastq-scan, ncbi-genome-download, pigz
 *
 * @modules gather, csvtk_concat
 *
 * @input record(meta, r1_files, r2_files, se_files, lr_files)
 * - `meta`: Groovy Map containing sample information
 * - `r1_files`: Illumina R1 read files (Set, for merging multiple runs)
 * - `r2_files`: Illumina R2 read files (Set, for merging multiple runs)
 * - `se_files`: Single-end read files (Set, for merging multiple runs)
 * - `lr_files`: Long read files (ONT/PacBio) or assembly for simulation
 *
 * @output sample_outputs
 * - `tsv`: A tab-delimited metadata file describing the valid samples
 *
 * @output run_outputs
 * - `csv`: Aggregated metadata from all samples
 */
nextflow.preview.types = true

include { GATHER as GATHER_MODULE } from '../../../modules/bactopia/gather/main'
include { CSVTK_CONCAT            } from '../../../modules/csvtk/concat/main'
include { filterWithData          } from 'plugin/nf-bactopia'
include { gatherCsvtk                  } from 'plugin/nf-bactopia'

workflow GATHER {
    take:
    samples: Channel<Record>

    main:
    GATHER_MODULE(samples)
    CSVTK_CONCAT(gatherCsvtk(GATHER_MODULE.out, 'tsv', [name: 'meta']), 'tsv', 'tsv')

    emit:
    // Downstream inputs
    reads = filterWithData(GATHER_MODULE.out, ['r1', 'r2', 'se', 'lr', 'fna'])

    // Published outputs
    sample_outputs = GATHER_MODULE.out
    run_outputs = CSVTK_CONCAT.out
}
