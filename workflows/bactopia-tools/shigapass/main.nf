#!/usr/bin/env nextflow
/**
 * Prediction of Shigella serotypes and differentiation from EIEC.
 *
 * This Bactopia Tool uses [ShigaPass](https://github.com/imanyass/ShigaPass) for in silico
 * prediction of serotypes in *Shigella* assemblies and differentiation between *Shigella*,
 * EIEC (Enteroinvasive *E. coli*) and non-*Shigella EIEC strains. ShigaPass analyzes
 * key antigenic determinants including O-antigen processing genes and invasion plasmid
 * antigens (ipa genes) to provide accurate serotype predictions following the
 * White-Kauffmann-Le Minor scheme. This enables rapid serological characterization
 * essential for epidemiological investigations and outbreak tracking.
 *
 * @status stable
 * @keywords Shigella, EIEC, serotyping, O-antigen, ipa genes, epidemiology, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, shigapass
 *
 * @subworkflows bactopiatool_init, shigapass
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.csv                           CSV file containing predicted Shigella or EIEC serotype for each sample
 *
 * @section Merged Results
 * @publish shigapass.csv                   Merged CSV file containing ShigaPass results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                          Tool execution logs including ShigaPass output
 * @publish logs/nf-*                        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                     Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { SHIGAPASS           } from '../../../subworkflows/shigapass/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SHIGAPASS(BACTOPIATOOL_INIT.out.assembly)

    publish:
    // Per-sample
    sample_outputs = SHIGAPASS.out.sample_outputs
    sample_nf_logs = collectNextflowLogs(SHIGAPASS.out.sample_outputs)
    // Run-level
    run_outputs = SHIGAPASS.out.run_outputs
    run_nf_logs = collectNextflowLogs(SHIGAPASS.out.run_outputs)
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
