#!/usr/bin/env nextflow
/**
 * Penicillin Binding Protein (PBP) typing for Streptococcus pneumoniae.
 *
 * This Bactopia Tool uses [pbptyper](https://github.com/rpetit3/pbptyper) for typing
 * the Penicillin Binding Protein (PBP) of _Streptococcus pneumoniae_ assemblies
 * to predict penicillin susceptibility.
 *
 * @status stable
 * @keywords streptococcus pneumoniae, pbp, penicillin resistance, typing, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing,resistance
 * @citation pbptyper
 *
 * @subworkflows bactopiatool_init, pbptyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            PBP typing results for each sample
 *
 * @section Merged Results
 * @publish pbptyper.tsv      Merged TSV file containing PBP typing results from all samples
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml      Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { PBPTYPER            } from '../../../subworkflows/pbptyper/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    PBPTYPER(BACTOPIATOOL_INIT.out.assembly)

    publish:
    // Per-sample
    sample_outputs = PBPTYPER.out.sample_outputs
    sample_nf_logs = collectNextflowLogs(PBPTYPER.out.sample_outputs)
    // Run-level
    run_outputs = PBPTYPER.out.run_outputs
    run_nf_logs = collectNextflowLogs(PBPTYPER.out.run_outputs)
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
