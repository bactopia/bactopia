#!/usr/bin/env nextflow
nextflow.preview.output = true

include { BACTOPIATOOL_INIT } from '../../lib/nf/bactopia-tool/init'
include { SISTR } from '../../subworkflows/sistr/main'

workflow {
    // Initialize samples
    samples = BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    SISTR(samples)

    publish {
        SISTR.out.tsv >> 'sistr'
        SISTR.out.merged_tsv >> 'sistr'
        SISTR.out.allele_fasta >> 'sistr'
        SISTR.out.allele_json >> 'sistr'
        SISTR.out.cgmlst_csv >> 'sistr'
        SISTR.out.logs >> 'logs/sistr'
        SISTR.out.nf_logs >> 'logs/nextflow'
    }

    output {
        meta {
            output_dir = "bactopia-tools/${meta.id}/${meta.output_dir}"
            logs_dir = "bactopia-tools/${meta.id}/${meta.logs_dir}"
            process_name = "${meta.process_name}"
        }
    }
}
