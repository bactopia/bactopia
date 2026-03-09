#!/usr/bin/env nextflow
/**
 * Search against translated nucleotide databases using translated nucleotide queries.
 *
 * This Bactopia Tool uses [TBLASTX](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs)
 * to query translated nucleotide sequences against translated nucleotide databases for
 * comprehensive homology search at the protein level.
 *
 * @status stable
 * @keywords fasta, blast, alignment, protein, translation, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,alignment,similarity-search
 * @citation blast
 *
 * @subworkflows bactopiatool_init, tblastx
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input tblastx_query
 * Path to nucleotide query sequence in FASTA format
 *
 * @input tblastx_db
 * Path to TBLASTX database for searching
 *
 * @section Per-Sample Results
 * @publish *.tblastx.tsv      TBLASTX alignment results in tabular format
 * @publish *.tblastx.html     Interactive HTML report of TBLASTX results
 *
 * @section Merged Results
 * @publish tblastx.tsv        Merged TBLASTX results from all samples
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
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    tblastx_query : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { TBLASTX           } from '../../../subworkflows/tblastx/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    TBLASTX(BACTOPIATOOL_INIT.out.assembly, params.tblastx_query)
    ch_sample_nf_logs = TBLASTX.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = TBLASTX.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    publish:
    sample_outputs = TBLASTX.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = TBLASTX.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}
output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
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

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
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
