#!/usr/bin/env nextflow
/**
 * Search against protein BLAST databases using protein queries.
 *
 * This Bactopia Tool uses [BLASTP](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs)
 * to query protein sequences against protein databases for sequence similarity search.
 * BLASTP compares a protein query to a protein database to find similar sequences.
 *
 * @status stable
 * @keywords fasta, blast, alignment, protein, similarity, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,alignment,similarity-search
 * @citation blast
 *
 * @subworkflows bactopiatool_init, blastp
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input blastp_query
 * Path to protein query sequence in FASTA format
 *
 * @input blastp_db
 * Path to BLAST protein database for searching
 *
 * @section Per-Sample Results
 * @publish *.blastp.tsv        BLASTP alignment results in tabular format
 * @publish *.blastp.html       Interactive HTML report of BLASTP results
 *
 * @section Merged Results
 * @publish blastp.tsv          Merged BLASTP results from all samples
 *
 * @section Execution Logs
 * @publish logs/**            Tool execution logs
 * @publish logs/nf-*          Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml        Software version information
   */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    blastp_query : Path
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { BLASTP              } from '../../../subworkflows/blastp/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    BLASTP(BACTOPIATOOL_INIT.out.blastdb, params.blastp_query)

    publish:
    // Per-sample
    sample_outputs = BLASTP.out.sample_outputs
    sample_nf_logs = collectNextflowLogs(BLASTP.out.sample_outputs)
    // Run-level
    run_outputs = BLASTP.out.run_outputs
    run_nf_logs = collectNextflowLogs(BLASTP.out.run_outputs)
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
