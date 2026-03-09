#!/usr/bin/env nextflow
/**
 * Sequence variation analysis of mcr-1 genes (mobilized colistin resistance).
 *
 * This Bactopia Tool uses [mcroni](https://github.com/liampshaw/mcroni) to identify _mcr-1_ genes in
 * assemblies and report sequence variations. If _mcr-1_ is found, the variations will be reported
 * and available in an output FASTA file.
 *
 * @status stable
 * @keywords mcr-1, colistin resistance, antimicrobial resistance, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,resistance
 * @citation mcroni
 *
 * @subworkflows bactopiatool_init, mcroni
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt            Summary of mcr-1 variants found
 * @publish *.fasta          FASTA file of mcr-1 variants
 *
 * @section Merged Results
 * @publish mcroni.tsv        Merged TSV file containing mcroni results from all samples
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
include { MCRONI            } from '../../../subworkflows/mcroni/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    MCRONI(BACTOPIATOOL_INIT.out.assembly)
    ch_sample_nf_logs = MCRONI.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MCRONI.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    publish:
    sample_outputs = MCRONI.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = MCRONI.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}
output {
    sample_outputs {
        path { r ->
            r.results  >> "${r.meta.output_dir}/"
            r.logs     >> "${r.meta.logs_dir}/"
            r.versions >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
    run_outputs {
        path { r ->
            r.results  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
