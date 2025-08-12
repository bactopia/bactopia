#!/usr/bin/env nextflow
nextflow.preview.output = true

include { BACTOPIATOOL_INIT } from '../../lib/nf/bactopia-tool/init'
include { SPATYPER } from '../../subworkflows/spatyper/main'

workflow {
    // Initialize samples
    samples = BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    SPATYPER(samples)

    publish {
        SPATYPER.out.tsv >> 'spatyper'
        SPATYPER.out.merged_tsv >> 'spatyper'
        SPATYPER.out.logs >> 'logs/spatyper'
        SPATYPER.out.nf_logs >> 'logs/nextflow'
    }

    output {
        meta {
            output_dir = "bactopia-tools/${meta.id}/${meta.output_dir}"
            logs_dir = "bactopia-tools/${meta.id}/${meta.logs_dir}"
            process_name = "${meta.process_name}"
        }
    }
}
