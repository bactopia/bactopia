#!/usr/bin/env nextflow
/**
 * Identify cap locus serotype and structure in Haemophilus influenzae assemblies.
 *
 * This Bactopia Tool uses [hicap](https://github.com/scwatts/hicap) with assemblies for
 * _in silico_ typing of the _Haemophilus influenzae_ capsular locus.
 *
 * @status stable
 * @keywords haemophilus influenzae, serotyping, capsular locus, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation hicap
 *
 * @subworkflows utils_bactopia-tools, hicap
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input hicap_database_dir
 * Path to hicap database directory
 *
 * @input hicap_model_fp
 * Path to hicap model file
 *
 * @section Per-Sample Results
 * @publish *.summary       Summary of serotype prediction
 * @publish *.gff           Annotated capsular locus in GFF format
 *
 * @section Merged Results
 * @publish hicap.tsv        Merged TSV file containing hicap results from all samples
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

    // Tool-specific parameters
    hicap_database_dir : Path?
    hicap_model_fp     : Path?
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { HICAP               } from '../../../subworkflows/hicap/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_hicap = HICAP(
        ch_bactopiatool.assembly,
        params.hicap_database_dir,
        params.hicap_model_fp
    )

    publish:
    // Per-sample
    sample_outputs = ch_hicap.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_hicap.sample_outputs)
    // Run-level
    run_outputs = ch_hicap.run_outputs
    run_nf_logs = collectNextflowLogs(ch_hicap.run_outputs)
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
