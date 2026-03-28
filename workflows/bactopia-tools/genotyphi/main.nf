#!/usr/bin/env nextflow
/**
 * Salmonella Typhi genotyping with lineage assignment.
 *
 * This Bactopia Tool uses [GenoTyphi](https://github.com/typhoidgenomics/genotyphi) to
 * call Typhi lineages, AMR determinants, and plasmid markers in Salmonella Typhi samples.
 * Samples are first processed by [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) using `mykrobe predict`
 * with `typhi` specified as the species. Then the Mykrobe results are processed by the
 * GenoTyphi parse_typhi_mykrobe.py script.
 *
 * @status stable
 * @keywords salmonella typhi, genotyping, lineage, amr, mykrobe, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,genotyping,lineage
 * @citation genotyphi, mykrobe
 *
 * @subworkflows bactopiatool_init, genotyphi
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.json           Mykrobe prediction results
 * @publish *.tsv            GenoTyphi parsed results
 *
 * @section Merged Results
 * @publish genotyphi.tsv     Merged TSV file containing GenoTyphi results from all samples
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

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { GENOTYPHI         } from '../../../subworkflows/genotyphi/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    GENOTYPHI(BACTOPIATOOL_INIT.out.reads)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = GENOTYPHI.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = GENOTYPHI.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = GENOTYPHI.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = GENOTYPHI.out.run_outputs
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
