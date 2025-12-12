/**
 * Remove human contamination from sequencing reads for SRA submission.
 *
 * This subworkflow uses the [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to identify
 * and remove human reads from sequencing data. It first initializes a human reference database
 * and then scrubs the input reads to ensure they meet SRA submission requirements.
 *
 * @status stable
 * @keywords contamination, human, scrub, sra, sequencing, fastq
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation srahumanscrubber
 *
 * @modules srahumanscrubber_initdb, srahumanscrubber_scrub
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end sequencing reads from clinical samples that may contain human contamination
 *
 * @output scrubbed             FASTQ files with human reads removed (cleaned reads ready for SRA submission)
 * @output scrubbed_extra       Additional FASTQ files containing reads removed during scrubbing
 * @output scrub_report         Summary report of human read detection and removal statistics
 * @output scrub_special_report Detailed report including read classification and contamination metrics
 * @output results              Aggregated results channel containing all output files
 * @output logs                 Aggregated logs channel containing all execution logs
 * @output nf_logs              Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions             Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { SRAHUMANSCRUBBER_INITDB } from '../../modules/srahumanscrubber/initdb/main'
include { SRAHUMANSCRUBBER_SCRUB  } from '../../modules/srahumanscrubber/scrub/main'
include { flattenPaths            } from 'plugin/nf-bactopia'
include { gather                  } from 'plugin/nf-bactopia'

workflow SRAHUMANSCRUBBER {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>

    main:
    SRAHUMANSCRUBBER_INITDB()
    SRAHUMANSCRUBBER_SCRUB(reads, SRAHUMANSCRUBBER_INITDB.out.db)

    emit:
    // Individual outputs
    scrubbed: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrubbed
    scrubbed_extra: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrubbed_extra
    scrub_report: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    scrub_special_report: Channel<Tuple<Map, Path>> = SRAHUMANSCRUBBER_SCRUB.out.scrub_special_report

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        SRAHUMANSCRUBBER_SCRUB.out.scrubbed,
        SRAHUMANSCRUBBER_SCRUB.out.scrub_report
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([SRAHUMANSCRUBBER_SCRUB.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([SRAHUMANSCRUBBER_SCRUB.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([SRAHUMANSCRUBBER_SCRUB.out.versions])
}
