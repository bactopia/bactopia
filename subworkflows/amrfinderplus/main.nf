/**
 * Find antimicrobial resistance genes and point mutations.
 *
 * This subworkflow uses [AMRFinderPlus](https://github.com/ncbi/amr) to identify acquired antimicrobial
 * resistance genes and some point mutations in protein or assembled nucleotide sequences.
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial resistance, gene prediction
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent, aggregation
 * @citation amrfinderplus
 *
 * @modules csvtk_concat, amrfinderplus_run
 *
 * @input tuple(meta, nucleotide, protein, gff)
 * - `meta`: Groovy Map containing sample information
 * - `nucleotide`: Assembled contigs in FASTA format for AMR gene detection
 * - `protein`: Optional protein sequences in FASTA format (Path?)
 * - `gff`: Optional GFF3 annotation file (Path?)
 *
 * @input db
 * Path to the AMRFinderPlus database directory containing reference data for AMR gene detection.
 *
 * @output report          AMR gene detection results for each sample in TSV format
 * @output merged_tsv      Combined AMR detection results from all samples in a single TSV file
 * @output mutation_report Point mutations associated with antimicrobial resistance
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { AMRFINDERPLUS_RUN } from '../../modules/amrfinderplus/run/main'
include { CSVTK_CONCAT      } from '../../modules/csvtk/concat/main'
include { flattenPaths      } from 'plugin/nf-bactopia'
include { gather            } from 'plugin/nf-bactopia'

workflow AMRFINDERPLUS {
    take:
    fasta: Channel<Tuple<Map, Path, Path?, Path?>>
    db: Path

    main:
    AMRFINDERPLUS_RUN(fasta, db)
    CSVTK_CONCAT(gather(AMRFINDERPLUS_RUN.out.report, 'amrfinderplus'), 'tsv', 'tsv')

    emit:
    // Individual outputs
    report: Channel<Tuple<Map, Set<Path>>> = AMRFINDERPLUS_RUN.out.report
    merged_tsv: Channel<Tuple<Map, Set<Path>>> = CSVTK_CONCAT.out.csv
    mutation_report: Channel<Tuple<Map, Set<Path>>> = AMRFINDERPLUS_RUN.out.mutation_report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.report,
        AMRFINDERPLUS_RUN.out.mutation_report,
        CSVTK_CONCAT.out.csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        AMRFINDERPLUS_RUN.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
