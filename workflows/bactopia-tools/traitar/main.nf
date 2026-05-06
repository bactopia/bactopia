#!/usr/bin/env nextflow
/**
 * Predict phenotypic traits from microbial genomes
 *
 * This Bactopia Tool uses [Traitar](https://github.com/nick-youngblut/traitar3/) to predict phenotypic traits from microbial genomes.
 *
 * @status stable
 * @keywords phenotype, traits, pfam, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, traitar
 *
 * @subworkflows utils_bactopia-tools, traitar
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input traitar_db
 * Path to a Pfam-A HMM file (optional)
 *
 * @input download_traitar
 * Boolean flag to trigger automatic database download
 *
 * @section Per-Sample Results
 * @publish *.majority.tsv              Majority-vote combined phenotype trait predictions
 * @publish *.single_votes.tsv          Single-votes combined phenotype trait predictions
 * @publish supplemental/*              Supplemental Traitar output files
 *
 * @section Merged Results
 * @publish traitar-majority.tsv      Merged majority-vote phenotype predictions from all samples
 * @publish traitar-single.tsv        Merged single-vote phenotype predictions from all samples
 *
 * @section Execution Logs
 * @publish logs/traitar/*                Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.enable.types = true

params {
    rundir : String

    // Tool-specific parameters
    traitar_db       : Path?
    download_traitar : Boolean
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { TRAITAR             } from '../../../subworkflows/traitar/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_traitar = TRAITAR(ch_bactopiatool.assembly, params.traitar_db, params.download_traitar)

    publish:
    // Per-sample
    sample_outputs = ch_traitar.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_traitar.sample_outputs)
    // Run-level
    run_outputs = ch_traitar.run_outputs
    run_nf_logs = collectNextflowLogs(ch_traitar.run_outputs)
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
