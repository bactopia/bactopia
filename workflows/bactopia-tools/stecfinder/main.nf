#!/usr/bin/env nextflow
/**
 * Serotype identification of Shiga toxin-producing E. coli.
 *
 * This Bactopia Tool uses [STECFinder](https://github.com/LanLab/STECFinder) to identify
 * the serotype of Shiga toxin-producing _E. coli_ (STEC) from sequencing data.
 * STECFinder determines the serotype as well as the O-antigen and H-antigens.
 *
 * @status stable
 * @keywords stec, serotype, e coli, shiga toxin, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,serotyping
 * @citation stecfinder
 *
 * @subworkflows utils_bactopia-tools, stecfinder
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            Serotype identification results
 *
 * @section Merged Results
 * @publish stecfinder.tsv   Merged TSV file containing STECFinder results from all samples
 *
 * @section Execution Logs
 * @publish logs/**          Tool execution logs
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { STECFINDER          } from '../../../subworkflows/stecfinder/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_stecfinder = STECFINDER(ch_bactopiatool.assembly_reads)

    publish:
    // Per-sample
    sample_outputs = ch_stecfinder.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_stecfinder.sample_outputs)
    // Run-level
    run_outputs = ch_stecfinder.run_outputs
    run_nf_logs = collectNextflowLogs(ch_stecfinder.run_outputs)
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
