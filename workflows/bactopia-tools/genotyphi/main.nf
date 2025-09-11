#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { GENOTYPHI         } from '../../../subworkflows/genotyphi/main'
include { paramsHelp        } from 'plugin/nf-bactopia'
include { workflowSummary   } from 'plugin/nf-bactopia'

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
    GENOTYPHI(BACTOPIATOOL_INIT.out.samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = GENOTYPHI.out.tsv.mix(
        GENOTYPHI.out.merged_tsv,
        GENOTYPHI.out.csv,
        GENOTYPHI.out.json
    )
    logs = GENOTYPHI.out.logs
    logs2 = GENOTYPHI.out.logs2
    nf_logs = GENOTYPHI.out.nf_logs
    nf_logs2 = GENOTYPHI.out.nf_logs2
    versions = GENOTYPHI.out.versions
    versions2 = GENOTYPHI.out.versions2
}

output {
    results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    logs {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    logs2 {
        path { meta, _file -> "${meta.logs_dir2}/" }
    }
    nf_logs {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    nf_logs2 {
        path { meta, file -> {
            file >> "${meta.logs_dir2}/nf${file.name}"
        } }
    }
    versions {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    versions2 {
        path { meta, _file -> "${meta.logs_dir2}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
