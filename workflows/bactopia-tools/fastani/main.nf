#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools'
include { FASTANI            } from '../../../subworkflows/fastani/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { paramsHelp         } from 'plugin/nf-bactopia'
include { workflowSummary    } from 'plugin/nf-bactopia'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    // Check if help is requested
    if (params.help || params.help_all) {
        log.info paramsHelp()
        exit 0
    }

    // Initialize and execute the workflow
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    // Reference if applicable
    ch_reference = Channel.empty()
    if (params.fastani_reference) {
        ch_reference << tuple([id:file(params.fastani_reference).getSimpleName()], file(params.fastani_reference))
    }

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        NCBIGENOMEDOWNLOAD(params.accessions ? file(params.accessions) : [])
        ch_reference = ch_reference.mix(NCBIGENOMEDOWNLOAD.out.bactopia_tools)
    }

    // Add query if pairwise
    ch_query = BACTOPIATOOL_INIT.out.samples
    if (params.fastani_pairwise) {
        ch_reference = ch_reference.mix(BACTOPIATOOL_INIT.out.samples)
        ch_query = ch_reference
    }

    // Run FastANI
    FASTANI(ch_query, ch_reference)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = FASTANI.out.results
    logs = FASTANI.out.logs
    nf_logs = FASTANI.out.nf_logs
    versions = FASTANI.out.versions
}

output {
    results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    logs {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    nf_logs {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    versions {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
