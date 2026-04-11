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
 * @subworkflows utils_bactopia-tools, mashtree, ncbigenomedownload
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

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { MASHTREE            } from '../../../subworkflows/mashtree/main'
include { NCBIGENOMEDOWNLOAD  } from '../../../subworkflows/ncbigenomedownload/main'
include { gather              } from 'plugin/nf-bactopia'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_samples = ch_bactopiatool.assembly

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        ch_ncbigenomedownload = NCBIGENOMEDOWNLOAD(params.accessions ? file(params.accessions) : [])
        ch_samples = ch_samples.mix(ch_ncbigenomedownload.assemblies)
    }
    ch_mashtree = MASHTREE(gather(ch_samples, 'fna', [name: 'mashtree']))

    publish:
    // Per-sample
    sample_outputs = ch_mashtree.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_mashtree.sample_outputs)
    // Run-level
    run_outputs = ch_mashtree.run_outputs
    run_nf_logs = collectNextflowLogs(ch_mashtree.run_outputs)
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
