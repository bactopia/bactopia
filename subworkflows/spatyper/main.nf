/**
 * Predict spa types of Staphylococcus aureus from genome assemblies.
 *
 * This subworkflow uses [spaTyper](https://github.com/HCGB-IGTP/spaTyper) to predict
 * the spa types of *Staphylococcus aureus* strains from assembled genomes based on
 * the polymorphic X region of the protein A gene (spa). It processes each sample
 * individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords staphylococcus aureus, spa typing, protein a, mrsa
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation spatyper
 *
 * @modules csvtk_concat, spatyper
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input repeats
 * Optional custom repeats database for spa typing
 *
 * @input repeat_order
 * Optional custom repeat order file for spa typing
 *
 * @output tsv         Per-sample TSV files containing spa typing results
 * @output merged_tsv  Consolidated TSV file containing spa typing from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SPATYPER as SPATYPER_MODULE } from '../../modules/spatyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow SPATYPER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    repeats: Path?
    repeat_order: Path?

    main:
    SPATYPER_MODULE(assembly, repeats, repeat_order)
    CSVTK_CONCAT(gather(SPATYPER_MODULE.out.tsv, 'spatyper'), 'tsv', 'tsv')

    emit:
    // Individual output
    tsv: Channel<Tuple<Map, Path>> = SPATYPER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate output
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        SPATYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
