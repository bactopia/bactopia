#!/usr/bin/env nextflow
/**
 * Search against nucleotide BLAST databases using nucleotide queries.
 *
 * This Bactopia Tool uses [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs)
 * to query nucleotide sequences against nucleotide databases for sequence similarity search.
 * BLASTN finds regions of local similarity between nucleotide sequences.
 *
 * @status stable
 * @keywords fasta, blast, alignment, nucleotide, similarity, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,alignment,similarity-search
 *
 * @subworkflows bactopiatool_init, blastn
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input blastn_query
 * Path to nucleotide query sequence in FASTA format
 *
 * @input blastn_db
 * Path to BLAST nucleotide database for searching
 *
 * @section Per-Sample Results
 * @publish *.blastn.tsv        BLASTN alignment results in tabular format
 * @publish *.blastn.html       Interactive HTML report of BLASTN results
 *
 * @section Merged Results
 * @publish blastn.tsv          Merged BLASTN results from all samples
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
    rundir   : String

    // Tool-specific parameters
    blastn_query : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { BLASTN            } from '../../../subworkflows/blastn/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    BLASTN(BACTOPIATOOL_INIT.out.assembly, params.blastn_query)
    ch_sample_nf_logs = BLASTN.out.sample_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    ch_run_nf_logs = BLASTN.out.run_outputs.flatMap { r -> r.nf_logs.collect { f -> tuple(r.meta, f) } }
    publish:
    sample_outputs = BLASTN.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = BLASTN.out.run_outputs
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
