#!/usr/bin/env nextflow
nextflow.preview.output = true

include { BACTOPIATOOL_INIT } from '../../lib/nf/bactopia-tool/init'
include { SHIGATYPER } from '../../subworkflows/shigatyper/main'

workflow {
    // Initialize samples
    samples = BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    SHIGATYPER(samples)

    publish {
        SHIGATYPER.out.tsv >> 'shigatyper'
        SHIGATYPER.out.merged_tsv >> 'shigatyper'
        SHIGATYPER.out.logs >> 'logs/shigatyper'
        SHIGATYPER.out.nf_logs >> 'logs/nextflow'
    }

    output {
        meta {
            output_dir = "bactopia-tools/${meta.id}/${meta.output_dir}"
            logs_dir = "bactopia-tools/${meta.id}/${meta.logs_dir}"
            process_name = "${meta.process_name}"
        }
    }
}
