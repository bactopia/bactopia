#!/usr/bin/env nextflow
/**
 * Search against protein BLAST databases using translated nucleotide queries.
 *
 * This Bactopia Tool uses [BLASTX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs)
 * to query translated nucleotide sequences against protein databases for protein homology search.
 *
 * @status stable
 * @keywords fasta, blast, alignment, protein, translation, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,alignment,similarity-search
 * @citation blast
 *
 * @subworkflows bactopiatool_init, blastx
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input blastx_query
 * Path to nucleotide query sequence in FASTA format
 *
 * @section Per-Sample Results
 * @publish *.blastx.tsv      BLASTX alignment results in tabular format
 * @publish *.blastx.html     Interactive HTML report of BLASTX results
 *
 * @section Merged Results
 * @publish blastx.tsv        Merged BLASTX results from all samples
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
    rundir   : String

    // Tool-specific parameters
    blastx_query : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { BLASTX            } from '../../../subworkflows/blastx/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    BLASTX(BACTOPIATOOL_INIT.out.blastdb, params.blastx_query)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = BLASTX.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = BLASTX.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = BLASTX.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = BLASTX.out.run_outputs
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
