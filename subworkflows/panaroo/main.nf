/**
 * Build a pangenome from GFF3 annotations using Panaroo.
 *
 * This subworkflow creates a pangenome from bacterial genome annotations using [Panaroo](https://github.com/gtonkinhill/panaroo).
 * Panaroo is a pangenome pipeline that produces polished pangenomes by removing errors and
 * contamination from input annotations. It generates gene presence/absence matrices and core-genome
 * alignments suitable for downstream phylogenetic analysis.
 *
 * @status stable
 * @keywords pangenome, pan-genome, comparative genomics, core-genome, alignment
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation panaroo
 *
 * @modules panaroo_run
 *
 * @input tuple(meta, gff)
 * - `meta`: Metadata map containing sample information including sample ID, name, and other attributes
 * - `gff`: Set of GFF3 annotation files representing the genomic annotations for each sample
 *
 * @output csv          Gene presence/absence matrix in CSV format showing which genes are present in each genome
 * @output aln          Core-genome alignment file containing genes present across all input genomes
 * @output filtered_aln Filtered core-genome alignment with recombination regions removed (if enabled)
 * @output panaroo_csv  Detailed Panaroo output including gene cluster information and statistics
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PANAROO_RUN  } from '../../modules/panaroo/run/main'
include { flattenPaths } from 'plugin/nf-bactopia'
include { gather       } from 'plugin/nf-bactopia'

workflow PANAROO {
    take:
    gff: Channel<Tuple<Map, Set<Path>>>

    main:
    PANAROO_RUN(gff)

    emit:
    // Individual outputs
    csv: Channel<Tuple<Map, Set<Path>>> = PANAROO_RUN.out.csv
    aln: Channel<Tuple<Map, Set<Path>>> = PANAROO_RUN.out.aln
    filtered_aln: Channel<Tuple<Map, Set<Path>>> = PANAROO_RUN.out.filtered_aln
    panaroo_csv: Channel<Tuple<Map, Set<Path>>> = PANAROO_RUN.out.panaroo_csv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PANAROO_RUN.out.supplemental,
        PANAROO_RUN.out.csv,
        PANAROO_RUN.out.aln,
        PANAROO_RUN.out.filtered_aln,
        PANAROO_RUN.out.panaroo_csv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([PANAROO_RUN.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([PANAROO_RUN.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([PANAROO_RUN.out.versions])
}
