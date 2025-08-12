#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { AMRFINDERPLUS     } from '../../../subworkflows/amrfinderplus/main'
include { DATASETS          } from '../../../modules/bactopia/datasets/main'
include { workflowSummary   } from 'plugin/nf-bactopia'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)
    
    DATASETS()
    AMRFINDERPLUS(BACTOPIATOOL_INIT.out.samples, DATASETS.out.amrfinderplus_db)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = AMRFINDERPLUS.out.report.mix(
        AMRFINDERPLUS.out.merged_tsv,
        AMRFINDERPLUS.out.mutation_report
    )
    logs = AMRFINDERPLUS.out.logs
    nf_logs = AMRFINDERPLUS.out.nf_logs
    versions = AMRFINDERPLUS.out.versions.mix(DATASETS.out.versions)
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
