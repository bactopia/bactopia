/**
 * Perform capsular typing of Streptococcus pneumoniae from NGS data.
 *
 * This subworkflow uses [PneumoCaT](https://github.com/ukhsa-collaboration/PneumoCaT) to
 * identify serotype-specific capsular loci and determine serotypes from next-generation
 * sequencing data. It provides comprehensive serotype determination including coverage
 * statistics and confidence scores for each sample.
 *
 * @status stable
 * @keywords streptococcus pneumoniae, serotype, capsular typing, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent
 * @citation pneumocat
 *
 * @modules pneumocat
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end reads in FASTQ format
 *
 * @output xml          Per-sample PneumoCaT detailed results in XML format with coverage information
 * @output txt          Per-sample summary reports with serotype calls and statistics
 * @output results      Aggregated results channel containing all output files
 * @output logs         Aggregated logs channel containing all execution logs
 * @output nf_logs      Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions     Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { PNEUMOCAT as PNEUMOCAT_MODULE } from '../../modules/pneumocat/main'
include { flattenPaths                  } from 'plugin/nf-bactopia'
include { gather                        } from 'plugin/nf-bactopia'

workflow PNEUMOCAT {
    take:
    reads: Channel<Tuple<Map, Set<Path>>>

    main:
    PNEUMOCAT_MODULE(reads)

    emit:
    // Individual outputs
    xml: Channel<Tuple<Map, Set<Path>>> = PNEUMOCAT_MODULE.out.xml
    txt: Channel<Tuple<Map, Set<Path>>> = PNEUMOCAT_MODULE.out.txt

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([
        PNEUMOCAT_MODULE.out.xml,
        PNEUMOCAT_MODULE.out.txt
    ])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.nf_logs])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([PNEUMOCAT_MODULE.out.versions])
}
