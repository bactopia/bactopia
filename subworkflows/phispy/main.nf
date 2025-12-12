/**
 * Prediction of prophages from bacterial genomes.
 *
 * This subworkflow identifies prophages in bacterial genomes using [PhiSpy](https://github.com/linsalrob/PhiSpy),
 * which combines similarity-based and composition-based strategies for accurate detection.
 * The tool identifies integrated phage sequences, extracts bacterial and phage regions,
 * and provides comprehensive annotation including GFF format for downstream analysis.
 *
 * @status stable
 * @keywords prophage, phage, bacterial, genome, mobile genetic elements
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation phispy
 *
 * @modules csvtk_concat, phispy
 *
 * @input gbk
 * Annotated bacterial genomes in GenBank format for prophage prediction
 *
 * @output tsv             PhiSpy prophage prediction results in TSV format
 * @output merged_tsv      Combined TSV file containing prophage results from all samples
 * @output information     Detailed information about identified prophages
 * @output bacteria_fasta  Extracted bacterial genomic sequences in FASTA format
 * @output bacteria_gbk    Extracted bacterial genomic sequences in GenBank format
 * @output phage_fasta     Extracted prophage sequences in FASTA format
 * @output phage_gbk       Extracted prophage sequences in GenBank format
 * @output prophage_gff    Prophage annotations in GFF3 format
 * @output prophage_tbl    Prophage protein table file
 * @output prophage_tsv    Prophage detailed annotations in TSV format
 * @output results         Aggregated results channel containing all output files
 * @output logs            Aggregated logs channel containing all execution logs
 * @output nf_logs         Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions        Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PHISPY as PHISPY_MODULE } from '../../modules/phispy/main'
include { CSVTK_CONCAT            } from '../../modules/csvtk/concat/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow PHISPY {
    take:
    gbk: Channel<Tuple<Map, Set<Path>>>

    main:
    PHISPY_MODULE(gbk)
    CSVTK_CONCAT(gather(PHISPY_MODULE.out.tsv, 'phispy'), 'tsv', 'tsv')

    emit:
    tsv: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.tsv
    merged_tsv: Channel<Tuple<Map, Path>> = CSVTK_CONCAT.out.csv
    information: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.information
    bacteria_fasta: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.bacteria_fasta
    bacteria_gbk: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.bacteria_gbk
    phage_fasta: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.phage_fasta
    phage_gbk: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.phage_gbk
    prophage_gff: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.prophage_gff
    prophage_tbl: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.prophage_tbl
    prophage_tsv: Channel<Tuple<Map, Path>> = PHISPY_MODULE.out.prophage_tsv

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.tsv,
        CSVTK_CONCAT.out.csv,
        PHISPY_MODULE.out.information,
        PHISPY_MODULE.out.bacteria_fasta,
        PHISPY_MODULE.out.bacteria_gbk,
        PHISPY_MODULE.out.phage_fasta,
        PHISPY_MODULE.out.phage_gbk,
        PHISPY_MODULE.out.prophage_gff,
        PHISPY_MODULE.out.prophage_tbl,
        PHISPY_MODULE.out.prophage_tsv
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.logs,
        CSVTK_CONCAT.out.logs
    ])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.nf_logs,
        CSVTK_CONCAT.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([
        PHISPY_MODULE.out.versions,
        CSVTK_CONCAT.out.versions
    ])
}
