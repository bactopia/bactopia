#!/usr/bin/env nextflow
nextflow.preview.output = true

include { BACTOPIATOOL_INIT } from '../../lib/nf/bactopia-tool/init'
include { SNIPPY } from '../../subworkflows/snippy/main'

workflow {
    // Initialize samples
    samples = BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    // Create reference channel
    ch_reference = Channel.fromPath(params.reference)
                          .map { file -> [[id: file.getSimpleName()], file] }

    SNIPPY(samples, ch_reference)

    publish {
        SNIPPY.out.consensus_subs_fa >> 'snippy'
        SNIPPY.out.consensus_subs_masked_fa >> 'snippy'
        SNIPPY.out.coverage >> 'snippy'
        SNIPPY.out.csv >> 'snippy'
        SNIPPY.out.filt_vcf >> 'snippy'
        SNIPPY.out.gff >> 'snippy'
        SNIPPY.out.html >> 'snippy'
        SNIPPY.out.raw_vcf >> 'snippy'
        SNIPPY.out.subs_vcf >> 'snippy'
        SNIPPY.out.tab >> 'snippy'
        SNIPPY.out.txt >> 'snippy'
        SNIPPY.out.vcf >> 'snippy'
        SNIPPY.out.logs >> 'logs/snippy'
        SNIPPY.out.nf_logs >> 'logs/nextflow'
    }

    output {
        meta {
            output_dir = "bactopia-tools/${meta.id}/${meta.output_dir}"
            logs_dir = "bactopia-tools/${meta.id}/${meta.logs_dir}"
            process_name = "${meta.process_name}"
        }
    }
}
