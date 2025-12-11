#!/usr/bin/env nextflow
nextflow.preview.types = true
/**
 * Bactopia Tool: Fastani.
 *
 * fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
 * The `fastani` module uses [FastANI](https://github.com/ParBLiSS/FastANI) to calculate the average
 * nucleotide identity (ANI) between your samples.
 * Although, sometimes you might be more interested in calculating the ANI of your samples against
 * a reference genome. Fortunately, using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download),
 * the `fastani` module allows you specify either a specific NCBI Assembly RefSeq accession (`--accession`)
 * or a species name (`--species`) for which to download all RefSeq genomes.
 *
 * @status stable
 *
 * @subworkflows bactopiatool_init, fastani
 *
 * @input rundir
 * Run directory containing Bactopia results
 *
 * @section Per-Sample Results
 * @publish *    Analysis results
 *
 * @section Merged Results
 * @publish merged-*    Aggregated results from all samples
 *
 * @section Execution Logs
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution logs
 *
 * @section Versions
 * @publish versions.yml Software version information
   */

params {
    rundir   : String

    // Tool-specific parameters
    fastani_reference : Path?
    fastani_pairwise  : Boolean
    species           : String?
    accession         : String?
    accessions        : Path?
}

include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools/main'
include { FASTANI            } from '../../../subworkflows/fastani/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()

    // Reference if applicable
    ch_reference = channel.empty() as Channel<Tuple<Map, Path>>
    if (params.fastani_reference) {
        ch_reference.mix(
            channel.of(tuple([id:params.fastani_reference.getSimpleName()], params.fastani_reference))
        )
    }

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        NCBIGENOMEDOWNLOAD(params.accessions)
        ch_reference = ch_reference.mix(NCBIGENOMEDOWNLOAD.out.bactopia_tools)
    }

    // Add query if pairwise
    ch_query = BACTOPIATOOL_INIT.out.samples
    if (params.fastani_pairwise) {
        ch_reference = ch_reference.mix(ch_query)
        ch_query = ch_reference
    }

    // Run FastANI
    FASTANI(ch_query, ch_reference)

    // Collect outputs
    ch_results = ch_results.mix(FASTANI.out.results)
    ch_logs = ch_logs.mix(FASTANI.out.logs)
    ch_nf_logs = ch_nf_logs.mix(FASTANI.out.nf_logs)
    ch_versions = ch_versions.mix(FASTANI.out.versions)

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
