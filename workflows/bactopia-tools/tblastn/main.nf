#!/usr/bin/env nextflow
/**
 * Search against translated nucleotide databases using protein queries.
 *
 * This Bactopia Tool uses [TBLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs)
 * to query protein sequences against translated nucleotide databases (contigs) for homology search.
 *
 * @status stable
 * @keywords fasta, blast, alignment, protein, nucleotide, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,alignment,similarity-search
 * @citation blast
 *
 * @subworkflows bactopiatool_init, tblastn
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input tblastn_query
 * Path to protein query sequence in FASTA format
 *
 * @input tblastn_db
 * Path to TBLASTN database for searching
 *
 * @section Per-Sample Results
 * @publish *.tblastn.tsv      TBLASTN alignment results in tabular format
 * @publish *.tblastn.html     Interactive HTML report of TBLASTN results
 *
 * @section Merged Results
 * @publish tblastn.tsv        Merged TBLASTN results from all samples
 *
 * @section Execution Logs
 * @publish logs/**            Tool execution logs
 * @publish logs/nf-*          Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml       Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    tblastn_query : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { TBLASTN           } from '../../../subworkflows/tblastn/main'

include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    TBLASTN(BACTOPIATOOL_INIT.out.blastdb, params.tblastn_query)

    ch_sample_nf_logs = collectNextflowLogs(TBLASTN.out.sample_outputs)
    ch_run_nf_logs = collectNextflowLogs(TBLASTN.out.run_outputs)

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = TBLASTN.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = TBLASTN.out.run_outputs
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
