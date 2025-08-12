#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { GENOTYPHI         } from '../../../subworkflows/genotyphi/main'
include { workflowSummary   } from 'plugin/nf-bactopia'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)
    
    GENOTYPHI(BACTOPIATOOL_INIT.out.samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    tsv = GENOTYPHI.out.tsv
    merged_tsv = GENOTYPHI.out.merged_tsv
    csv = GENOTYPHI.out.csv
    json = GENOTYPHI.out.json
    logs = GENOTYPHI.out.logs
    nf_logs = GENOTYPHI.out.nf_logs
    versions = GENOTYPHI.out.versions
}

output {
    tsv {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    merged_tsv {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    csv {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    json {
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
