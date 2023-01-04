#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
include { create_input_channel; check_input_fofn } from '../lib/nf/bactopia'
include { get_resources; get_schemas; print_efficiency } from '../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
ADAPTERS = params.adapters ? file(params.adapters) : file(params.empty_adapters)
PHIX = params.phix ? file(params.phix) : file(params.empty_phix)
SCHEMAS = get_schemas()

WorkflowMain.initialise(workflow, params, log, schema_filename=SCHEMAS)
run_type = WorkflowBactopia.initialise(workflow, params, log, schema_filename=SCHEMAS)
if (params.check_samples) {
    check_input_fofn()
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Core
include { GATHER_SAMPLES } from '../modules/local/bactopia/gather_samples/main'
include { QC_READS } from '../modules/local/bactopia/qc_reads/main'
include { SCRUBBER } from '../subworkflows/local/scrubber/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'  addParams( options: [publish_to_base: true] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow CLEANYERREADS {
    print_efficiency(RESOURCES.MAX_CPUS)
    ch_versions = Channel.empty()
    
    // Core steps
    GATHER_SAMPLES(create_input_channel(run_type, params.genome_size))
    ch_versions = ch_versions.mix(GATHER_SAMPLES.out.versions.first())

    if (params.enable_scrubber) {
        // Remove host reads
        SCRUBBER(GATHER_SAMPLES.out.raw_fastq)
        ch_versions = ch_versions.mix(SCRUBBER.out.versions)

        // Clean up scrubbed reads
        QC_READS(SCRUBBER.out.scrubbed.combine(Channel.fromPath(ADAPTERS)).combine(Channel.fromPath(PHIX)))
        ch_versions = ch_versions.mix(QC_READS.out.versions.first())
    } else {
        // Clean up raw reads
        QC_READS(GATHER_SAMPLES.out.raw_fastq.combine(Channel.fromPath(ADAPTERS)).combine(Channel.fromPath(PHIX)))
        ch_versions = ch_versions.mix(QC_READS.out.versions.first())
    }

    // Collect Versions
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile())
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/
workflow.onComplete {
    workDir = new File("${workflow.workDir}")

    println """
    Bactopia Execution Summary
    ---------------------------
    Bactopia Version : ${workflow.manifest.version}
    Nextflow Version : ${nextflow.version}
    Command Line     : ${workflow.commandLine}
    Resumed          : ${workflow.resume}
    Completed At     : ${workflow.complete}
    Duration         : ${workflow.duration}
    Success          : ${workflow.success}
    Exit Code        : ${workflow.exitStatus}
    Error Report     : ${workflow.errorReport ?: '-'}
    Launch Dir       : ${workflow.launchDir}
    """
}

/*
========================================================================================
    THE END
========================================================================================
*/
