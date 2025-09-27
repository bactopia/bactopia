#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools'
include { MASHTREE           } from '../../../subworkflows/mashtree/main'
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
    ch_samples = BACTOPIATOOL_INIT.out.samples

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        NCBIGENOMEDOWNLOAD(params.accessions ? file(params.accessions) : [])
        ch_samples = ch_samples.mix(NCBIGENOMEDOWNLOAD.out.bactopia_tools)
    }

    MASHTREE(ch_samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = MASHTREE.out.tree.mix(
        MASHTREE.out.matrix,
        MASHTREE.out.sketches
    )
    nf_logs = MASHTREE.out.nf_logs
    versions = MASHTREE.out.versions
}

output {
    results {
        path { meta, _file -> "${meta.output_dir}/" }
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
