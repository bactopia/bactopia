#!/usr/bin/env nextflow
nextflow.preview.output = true

include { BACTOPIATOOL_INIT } from '../../lib/nf/bactopia-tool/init'
include { TBPROFILER } from '../../subworkflows/tbprofiler/main'

workflow {
    // Initialize samples
    samples = BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    TBPROFILER(samples)

    publish {
        TBPROFILER.out.bam >> 'tbprofiler'
        TBPROFILER.out.csv >> 'tbprofiler'
        TBPROFILER.out.merged_csv >> 'tbprofiler'
        TBPROFILER.out.json >> 'tbprofiler'
        TBPROFILER.out.txt >> 'tbprofiler'
        TBPROFILER.out.vcf >> 'tbprofiler'
        TBPROFILER.out.logs >> 'logs/tbprofiler'
        TBPROFILER.out.nf_logs >> 'logs/nextflow'
    }

    output {
        meta {
            output_dir = "bactopia-tools/${meta.id}/${meta.output_dir}"
            logs_dir = "bactopia-tools/${meta.id}/${meta.logs_dir}"
            process_name = "${meta.process_name}"
        }
    }
}
