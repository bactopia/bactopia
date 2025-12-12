/**
 * Predict emm types of Streptococcus pyogenes from genome assemblies.
 *
 * This subworkflow uses [emmtyper](https://github.com/MDU-PHL/emmtyper) to predict
 * the emm types of *Streptococcus pyogenes* strains from assembled genomes. It processes
 * each sample individually and aggregates the results into a single consolidated report.
 *
 * @status stable
 * @keywords streptococcus pyogenes, emm typing, gas, m protein
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation, database-dependent
 * @citation emmtyper
 *
 * @modules emmtyper, csvtk_concat
 *
 * @input tuple(meta, assembly)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `assembly`: Set of assembled contigs in FASTA format to be analyzed for emm genes
 *
 * @input blastdb
 * Optional BLAST database containing emm gene reference sequences for improved typing accuracy
 *
 * @output tsv         Per-sample TSV files containing emm typing results
 * @output merged_tsv  Consolidated TSV file containing emm typing from all samples
 * @output results     Aggregated results channel containing all output files
 * @output logs        Aggregated logs channel containing all execution logs
 * @output nf_logs     Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions    Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { EMMTYPER as EMMTYPER_MODULE } from '../../modules/emmtyper/main'
include { CSVTK_CONCAT                } from '../../modules/csvtk/concat/main'
include { flattenPaths                } from 'plugin/nf-bactopia'
include { gather                      } from 'plugin/nf-bactopia'

workflow EMMTYPER {
    take:
    assembly: Channel<Tuple<Map, Set<Path>>>
    blastdb: Path?

    main:
    EMMTYPER_MODULE(assembly, blastdb)
    CSVTK_CONCAT(gather(EMMTYPER_MODULE.out.tsv, 'emmtyper'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    tsv: Channel<Tuple<Map, Path>> = EMMTYPER_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        EMMTYPER_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        EMMTYPER_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        CSVTK_CONCAT.out.nf_logs,
        EMMTYPER_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        EMMTYPER_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
