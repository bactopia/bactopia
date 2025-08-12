#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT    } from '../../../subworkflows/utils/bactopia-tools'
include { NCBIGENOMEDOWNLOAD   } from '../../../subworkflows/ncbigenomedownload/main'
include { workflowSummary      } from 'plugin/nf-bactopia'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    // Note: ncbigenomedownload doesn't need sample input, it downloads based on params
    NCBIGENOMEDOWNLOAD()

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = NCBIGENOMEDOWNLOAD.out.gbk.mix(
        NCBIGENOMEDOWNLOAD.out.fna,
        NCBIGENOMEDOWNLOAD.out.rm,
        NCBIGENOMEDOWNLOAD.out.features,
        NCBIGENOMEDOWNLOAD.out.gff,
        NCBIGENOMEDOWNLOAD.out.faa,
        NCBIGENOMEDOWNLOAD.out.gpff,
        NCBIGENOMEDOWNLOAD.out.wgs_gbk,
        NCBIGENOMEDOWNLOAD.out.cds,
        NCBIGENOMEDOWNLOAD.out.rna,
        NCBIGENOMEDOWNLOAD.out.rna_fna,
        NCBIGENOMEDOWNLOAD.out.report,
        NCBIGENOMEDOWNLOAD.out.stats
    )
    logs = NCBIGENOMEDOWNLOAD.out.logs
    nf_logs = NCBIGENOMEDOWNLOAD.out.nf_logs
    versions = NCBIGENOMEDOWNLOAD.out.versions
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
