#!/usr/bin/env nextflow
/**
 * Comprehensive typing of Staphylococcus aureus genomes.
 *
 * This Bactopia Tool is a subworkflow that includes multiple tools specific for typing
 * _Staphylococcus aureus_ features. Currently `staphtyper` includes:
 * 1. [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) - agr locus type and operon variants
 * 2. [spaTyper](https://github.com/HCGB-IGTP/spaTyper) - spa type
 * 3. [sccmec](https://github.com/rpetit3/sccmec) - SCCmec type
 *
 * @status stable
 * @keywords staphylococcus aureus, agr, spa, sccmec, typing, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing,workflow
 * @citation agrvate, spatyper, sccmec
 *
 * @subworkflows bactopiatool_init, staphtyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Comprehensive Typing
 * @note Results from all included typing tools
 * @publish staphtyper.tsv    Merged summary containing agr, spa, and SCCmec typing results
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs from all typing tools
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml        Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { STAPHTYPER          } from '../../../subworkflows/staphtyper/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_staphtyper = STAPHTYPER(
        ch_bactopiatool.assembly,
        params.spatyper_repeats,
        params.spatyper_repeat_order
    )

    publish:
    // Per-sample
    sample_outputs = ch_staphtyper.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_staphtyper.sample_outputs)
    // Run-level
    run_outputs = ch_staphtyper.run_outputs
    run_nf_logs = collectNextflowLogs(ch_staphtyper.run_outputs)
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
