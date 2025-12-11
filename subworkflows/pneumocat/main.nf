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
 * @modules pneumocat as pneumocat_module
 *
 * @input fastq
 * Channel containing fastq data
 *
 * @output xml      Xml
 * @output txt      Txt
 * @output results  Aggregated results channel containing all output files
 * @output logs     Aggregated logs channel containing all execution logs
 * @output nf_logs  Aggregated Nextflow execution logs from all processes
 * @output versions Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PNEUMOCAT as PNEUMOCAT_MODULE } from '../../modules/pneumocat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow PNEUMOCAT {
    take:
    fastq: Channel<Tuple<Map, Set<Path>>>

    main:
    PNEUMOCAT_MODULE(fastq)

    emit:
    // Individual outputs
    xml: Channel<Tuple<Map, Path>> = PNEUMOCAT_MODULE.out.xml
    txt: Channel<Tuple<Map, Path>> = PNEUMOCAT_MODULE.out.txt

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PNEUMOCAT_MODULE.out.xml,
        PNEUMOCAT_MODULE.out.txt
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.versions])
}
