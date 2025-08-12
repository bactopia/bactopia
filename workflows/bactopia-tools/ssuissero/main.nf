#!/usr/bin/env nextflow
nextflow.preview.output = true

include { BACTOPIATOOL_INIT } from '../../lib/nf/bactopia-tool/init'
include { SSUISSERO } from '../../subworkflows/ssuissero/main'

workflow {
    // Initialize samples
    samples = BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    SSUISSERO(samples)

    publish {
        SSUISSERO.out.tsv >> 'ssuissero'
        SSUISSERO.out.merged_tsv >> 'ssuissero'
        SSUISSERO.out.logs >> 'logs/ssuissero'
        SSUISSERO.out.nf_logs >> 'logs/nextflow'
    }

    output {
        meta {
            output_dir = "bactopia-tools/${meta.id}/${meta.output_dir}"
            logs_dir = "bactopia-tools/${meta.id}/${meta.logs_dir}"
            process_name = "${meta.process_name}"
        }
    }
}
