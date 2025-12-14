#!/usr/bin/env nextflow
/**
 * Rapid phylogenetic tree construction using Mash distances.
 *
 * This Bactopia Tool uses [Mashtree](https://github.com/lskatz/mashtree) to create a phylogenetic tree
 * of samples using [Mash](https://github.com/marbl/Mash) distances. It can include reference
 * genomes from RefSeq by downloading them with NCBI genome download.
 *
 * @status stable
 * @keywords phylogeny, tree, mash, distance, comparative genomics, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,phylogeny,comparative
 * @citation mashtree
 *
 * @subworkflows bactopiatool_init, mashtree, ncbigenomedownload
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input species
 * Species name to download all RefSeq genomes for comparison
 *
 * @input accession
 * Specific NCBI Assembly RefSeq accession to download
 *
 * @input accessions
 * Path to file containing list of NCBI accessions to download
 *
 * @section Phylogenetic Analysis
 * @publish mashtree.dnd          Newick format tree file
 * @publish mashtree.tsv          Tab-delimited distance matrix
 *
 * @section Merged Results
 * @publish mashtree-summary.tsv  Merged summary of all Mashtree results
 *
 * @section Execution Logs
 * @publish logs/**               Tool execution logs
 * @publish logs/nf-*             Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml          Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools/main'
include { MASHTREE           } from '../../../subworkflows/mashtree/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { formatSamples      } from 'plugin/nf-bactopia'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    ch_samples = BACTOPIATOOL_INIT.out.assembly

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        NCBIGENOMEDOWNLOAD(params.accessions ? file(params.accessions) : [])
        ch_samples = ch_samples.mix(NCBIGENOMEDOWNLOAD.out.bactopia_tools)
    }

    MASHTREE(ch_samples)

    // Collect outputs
    ch_results = ch_results.mix(MASHTREE.out.results)
    ch_logs = ch_logs.mix(MASHTREE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MASHTREE.out.nf_logs)
    ch_versions = ch_versions.mix(MASHTREE.out.versions)    

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
