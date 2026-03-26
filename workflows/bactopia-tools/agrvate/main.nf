#!/usr/bin/env nextflow
/**
 * Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.
 *
 * This Bactopia Tool uses [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) to rapidly
 * identify the _agr_ locus type and detect agr operon variants in _Staphylococcus aureus_ assemblies.
 * The agr system is a key quorum-sensing regulator of virulence in S. aureus.
 *
 * @status stable
 * @keywords staphylococcus aureus, agr, typing, virulence, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation agrvate
 *
 * @subworkflows bactopiatool_init, agrvate
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                 Agr locus typing results in TSV format
 * @publish supplemental/*        Supplemental output files including detailed agr analysis
 *
 * @section Merged Results
 * @publish agrvate.tsv            Combined agr typing results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                Tool execution logs
 * @publish logs/nf-*              Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { AGRVATE           } from '../../../subworkflows/agrvate/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    AGRVATE(BACTOPIATOOL_INIT.out.assembly)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = AGRVATE.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = AGRVATE.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = AGRVATE.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = AGRVATE.out.run_outputs
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
