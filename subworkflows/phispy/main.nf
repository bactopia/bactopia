/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @modules phispy as phispy_module, csvtk_concat
 *
 * @input gbk
 * Channel containing gbk data
 *
 * @output tsv            Tsv
 * @output merged_tsv     Merged Tsv
 * @output information    Information
 * @output bacteria_fasta Bacteria Fasta
 * @output bacteria_gbk   Bacteria Gbk
 * @output phage_fasta    Phage Fasta
 * @output phage_gbk      Phage Gbk
 * @output prophage_gff   Prophage Gff
 * @output prophage_tbl   Prophage Tbl
 * @output prophage_tsv   Prophage Tsv
 * @output results        Aggregated results channel containing all output files
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution logs from all processes
 * @output versions       Aggregated version information from all executed tools
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
