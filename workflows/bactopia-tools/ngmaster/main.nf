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
 * @subworkflows bactopiatool_init, ngmaster
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
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { NGMASTER          } from '../../../subworkflows/ngmaster/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    NGMASTER(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = NGMASTER.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = NGMASTER.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    publish:
    sample_outputs = NGMASTER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = NGMASTER.out.run_outputs
    run_nf_logs = ch_run_nf_logs
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
