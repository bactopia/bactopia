#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools'
include { SNIPPY        } from '../../../subworkflows/snippy/main'
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
    SNIPPY(BACTOPIATOOL_INIT.out.samples)

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = BACTOPIATOOL_INIT.out.aligned_fa.mix(
        BACTOPIATOOL_INIT.out.annotated_vcf,
        BACTOPIATOOL_INIT.out.bam,
        BACTOPIATOOL_INIT.out.bai,
        BACTOPIATOOL_INIT.out.bed,
        BACTOPIATOOL_INIT.out.consensus_fa,
        BACTOPIATOOL_INIT.out.consensus_subs_fa,
        BACTOPIATOOL_INIT.out.consensus_subs_masked_fa,
        BACTOPIATOOL_INIT.out.coverage,
        BACTOPIATOOL_INIT.out.csv,
        BACTOPIATOOL_INIT.out.filt_vcf,
        BACTOPIATOOL_INIT.out.gff,
        BACTOPIATOOL_INIT.out.html,
        BACTOPIATOOL_INIT.out.raw_vcf,
        BACTOPIATOOL_INIT.out.subs_vcf,
        BACTOPIATOOL_INIT.out.tab,
        BACTOPIATOOL_INIT.out.txt,
        BACTOPIATOOL_INIT.out.vcf
    )
    logs = SNIPPY.out.logs
    nf_logs = SNIPPY.out.nf_logs
    versions = SNIPPY.out.versions
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
