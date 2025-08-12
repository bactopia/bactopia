#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { MOBSUITE          } from '../../../subworkflows/mobsuite/main'
include { workflowSummary   } from 'plugin/nf-bactopia'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)
    
    MOBSUITE(BACTOPIATOOL_INIT.out.samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    chromosome = MOBSUITE.out.chromosome
    contig_report = MOBSUITE.out.contig_report
    plasmids = MOBSUITE.out.plasmids
    mobtyper_results = MOBSUITE.out.mobtyper_results
    merged_reports = MOBSUITE.out.merged_reports
    logs = MOBSUITE.out.logs
    nf_logs = MOBSUITE.out.nf_logs
    versions = MOBSUITE.out.versions
}

output {
    chromosome {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    contig_report {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    plasmids {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    mobtyper_results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    merged_reports {
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
