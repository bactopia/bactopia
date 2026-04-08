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
    rundir : String

    // Tool-specific parameters
    blastx_query : Value<Path>
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { BLASTX              } from '../../../subworkflows/blastx/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_blastx = BLASTX(ch_bactopiatool.blastdb, params.blastx_query)

    publish:
    // Per-sample
    sample_outputs = ch_blastx.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_blastx.sample_outputs)
    // Run-level
    run_outputs = ch_blastx.run_outputs
    run_nf_logs = collectNextflowLogs(ch_blastx.run_outputs)
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
