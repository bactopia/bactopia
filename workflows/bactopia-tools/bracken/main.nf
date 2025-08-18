#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { BRACKEN           } from '../../../subworkflows/bracken/main'
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
    if (params.help) {
        log.info paramsHelp()
        exit 0
    }

    // Initialize and execute the workflow
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)
    BRACKEN(BACTOPIATOOL_INIT.out.samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = BRACKEN.out.tsv.mix(
        BRACKEN.out.merged_tsv,
        BRACKEN.out.classified,
        BRACKEN.out.unclassified,
        BRACKEN.out.kraken2_report,
        BRACKEN.out.kraken2_output,
        BRACKEN.out.bracken_report,
        BRACKEN.out.abundances,
        BRACKEN.out.classification,
        BRACKEN.out.adjusted_abundances,
        BRACKEN.out.merged_adjusted_abundances
    )
    logs = BRACKEN.out.logs
    nf_logs = BRACKEN.out.nf_logs
    versions = BRACKEN.out.versions
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
