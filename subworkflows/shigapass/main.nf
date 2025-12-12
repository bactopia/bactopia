/**
 * Predict serotypes of Shigella from assemblies.
 *
 * This subworkflow uses [ShigaPass](https://github.com/imanyass/ShigaPass) to predict
 * serotypes of *Shigella* strains from assembled genomes. It analyzes the presence
 * and composition of antigen-encoding genes to classify isolates into their known serotypes.
 *
 * @status stable
 * @keywords shigella, serotype, typing, prediction, antigen genes
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation shigapass
 *
 * @modules shigapass, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv         Per-sample TSV files containing Shigella serotype predictions
 * @output merged_tsv  Consolidated TSV file containing serotype predictions from all samples
 * @output flex_tsv    Per-sample TSV files containing flexible serotype predictions
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SHIGAPASS as SHIGAPASS_MODULE } from '../../modules/shigapass/main'
include { CSVTK_CONCAT                  } from '../../modules/csvtk/concat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow SHIGAPASS {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>

    main:
    SHIGAPASS_MODULE(assembly)
    CSVTK_CONCAT(gather(SHIGAPASS_MODULE.out.tsv, 'shigapass'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = SHIGAPASS_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    flex_tsv: Channel<Tuple<Map, Path>> = SHIGAPASS_MODULE.out.flex_tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.tsv,
        SHIGAPASS_MODULE.out.flex_tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SHIGAPASS_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
