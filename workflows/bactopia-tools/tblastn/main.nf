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
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    tblastn_query : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { TBLASTN           } from '../../../subworkflows/tblastn/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    TBLASTN(
        BACTOPIATOOL_INIT.out.samples,
        params.tblastn_query
    )

    // Collect outputs
    ch_results = ch_results.mix(TBLASTN.out.results)
    ch_logs = ch_logs.mix(TBLASTN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(TBLASTN.out.nf_logs)
    ch_versions = ch_versions.mix(TBLASTN.out.versions)

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}
