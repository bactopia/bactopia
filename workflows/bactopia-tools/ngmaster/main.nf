#!/usr/bin/env nextflow
/**
 * Multi-antigen sequence typing of Neisseria gonorrhoeae.
 *
 * This Bactopia Tool uses [ngmaster](https://github.com/MDU-PHL/ngmaster) for
 * _in silico_ multi-antigen sequence typing (NG-MAST) of _Neisseria gonorrhoeae_.
 *
 * @status stable
 * @keywords neisseria gonorrhoeae, ng-mast, typing, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation ngmaster
 *
 * @subworkflows utils_bactopia-tools, ngmaster
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            NG-MAST typing results
 *
 * @section Merged Results
 * @publish ngmaster.tsv      Merged TSV file containing NG-MASTER results from all samples
 *
 * @section Execution Logs
 * @publish logs/**          Tool execution logs
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.enable.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { NGMASTER            } from '../../../subworkflows/ngmaster/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_ngmaster = NGMASTER(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_ngmaster.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_ngmaster.sample_outputs)
    // Run-level
    run_outputs = ch_ngmaster.run_outputs
    run_nf_logs = collectNextflowLogs(ch_ngmaster.run_outputs)
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
