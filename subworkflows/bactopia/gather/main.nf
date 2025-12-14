/**
 * Search, validate, gather, and standardize input samples.
 *
 * This subworkflow processes raw input samples through validation, standardization, and metadata
 * collection. It handles various input types including local FASTQ files, SRA/ENA accessions,
 * NCBI assembly accessions, and assemblies. The workflow can merge multiple sequencing runs,
 * download remote data, and simulate reads from assemblies using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art).
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1_files, r2_files, se_files, lr_files) with Set<Path> slots (pre-merge)
 * - Output: tuple(meta, r1, r2, se, lr) with Path? slots (post-merge, consolidated)
 *
 * @status stable
 * @keywords validation, download, merging, simulation, metadata, fastq, sra, ena, art
 * @tags complexity:complex input-type:single output-type:multiple features:aggregation, resource-download, conditional-logic
 * @citation art, fastq-dl, fastq-scan, ncbi-genome-download, pigz
 *
 * @modules gather, csvtk_concat
 *
 * @input tuple(meta, r1_files, r2_files, se_files, lr_files)
 * - `meta`: Groovy Map containing sample information
 * - `r1_files`: Illumina R1 read files (Set, for merging multiple runs)
 * - `r2_files`: Illumina R2 read files (Set, for merging multiple runs)
 * - `se_files`: Single-end read files (Set, for merging multiple runs)
 * - `lr_files`: Long read files (ONT/PacBio) or assembly for simulation
 *
 * @output tsv         Per-sample metadata files in TSV format
 * @output merged_tsv  Consolidated metadata file containing information from all samples
 * @output reads       Tuple with explicit read slots: (meta, r1, r2, se, lr) where each is Path?
 * @output error       Error messages from validation or download failures
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { GATHER as GATHER_MODULE } from '../../../modules/bactopia/gather/main'
include { CSVTK_CONCAT            } from '../../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow GATHER {
    take:
    samples: Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>, Set<Path>>>

    main:
    GATHER_MODULE(samples)

    // Merge meta values for each sample
    CSVTK_CONCAT(gather(GATHER_MODULE.out.tsv, 'meta'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = GATHER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>> = GATHER_MODULE.out.reads
    error: Channel<Tuple<Map, Set<Path>>> = GATHER_MODULE.out.error

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.tsv,
        GATHER_MODULE.out.error,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        GATHER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        GATHER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
