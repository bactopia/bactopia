/**
 * Predict Salmonella serotypes from genome assemblies.
 *
 * This subworkflow uses [SeqSero2](https://github.com/denglab/SeqSero2) to predict
 * the serotypes of *Salmonella* strains from assembled genomes. It processes each
 * sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords salmonella, serotype, prediction, foodborne, enteric
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation seqsero2
 *
 * @modules csvtk_concat, seqsero2
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing Salmonella serotype predictions
 * @output txt         Per-sample TXT files containing detailed serotype predictions
 * @output merged_tsv  Consolidated TSV file containing serotype predictions from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SEQSERO2 as SEQSERO2_MODULE } from '../../modules/seqsero2/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SEQSERO2 {
    take:
    seqs: Channel<Tuple<Map, Path>>

    main:
    SEQSERO2_MODULE(seqs)
    CSVTK_CONCAT(gather(SEQSERO2_MODULE.out.tsv, 'seqsero2'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Set<Path>>> = SEQSERO2_MODULE.out.tsv
    txt: Channel<Tuple<Map, Set<Path>>> = SEQSERO2_MODULE.out.txt
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.tsv,
        SEQSERO2_MODULE.out.txt,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SEQSERO2_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
