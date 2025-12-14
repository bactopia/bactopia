/**
 * Calculate Average Nucleotide Identity (ANI) between genomes.
 *
 * This subworkflow uses [FastANI](https://github.com/ParBLiSS/FastANI) to compute
 * whole-genome Average Nucleotide Identity (ANI) values between query genomes
 * and reference genomes. ANI is a robust measure of genomic similarity used for
 * species delineation in microbial taxonomy. The results are aggregated into
 * a single consolidated report.
 *
 * @status stable
 * @keywords ani, average nucleotide identity, taxonomy, species, comparison
 * @tags complexity:moderate input-type:multiple output-type:multiple features:aggregation
 * @citation fastani
 *
 * @modules csvtk_concat, fastani
 *
 * @input tuple(meta, query)
 * - `meta`: Groovy Map containing sample information
 * - `query`: Query genomes in FASTA format for ANI calculation
 *
 * @input tuple(meta, reference)
 * - `meta`: Groovy Map containing sample information
 * - `reference`: Reference genomes in FASTA format for ANI calculation
 *
 * @output tsv             Per-sample TSV files containing ANI values against reference genomes
 * @output merged_tsv      Consolidated TSV file containing ANI values from all comparisons
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { FASTANI as FASTANI_MODULE } from '../../modules/fastani/main'
include { CSVTK_CONCAT              } from '../../modules/csvtk/concat/main'
include { flattenPaths              } from 'plugin/nf-bactopia'
include { gather                    } from 'plugin/nf-bactopia'

workflow FASTANI {
    take:
    query: Channel<Tuple<Map, Path>>
    reference: Channel<Tuple<Map, Path>>

    main:
    FASTANI_MODULE(gather(query, 'query', 'fasta'), gather(reference, 'reference', 'fasta'))
    CSVTK_CONCAT(gather(FASTANI_MODULE.out.tsv, 'fastani'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = FASTANI_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        FASTANI_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
