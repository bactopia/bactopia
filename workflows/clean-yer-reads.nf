#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
include { create_input_channel; check_input_fofn } from '../lib/nf/bactopia'
include { get_schemas } from '../lib/nf/functions'

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
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
include { GATHER } from '../subworkflows/local/gather/main'
include { QC } from '../subworkflows/local/qc/main'
include { SCRUBBER } from '../subworkflows/local/scrubber/main'
include { CSVTK_CONCAT } from '../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'scrubber-concat', process_name: params.merge_folder] )
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
    ch_versions = Channel.empty()
    
    // Core steps
    GATHER(create_input_channel(run_type, params.genome_size, params.species))
    ch_versions = ch_versions.mix(GATHER.out.versions.first())

    if (params.use_k2scrubber || params.use_srascrubber) {
        // Remove host reads
        SCRUBBER(GATHER.out.fastq_only)
        ch_versions = ch_versions.mix(SCRUBBER.out.versions)

        // Merge scrub reports
        SCRUBBER.out.tsv.collect{meta, summary -> summary}.map{ summary -> [[id:'scrubber'], summary]}.set{ ch_merge_scrubber }
        CSVTK_CONCAT(ch_merge_scrubber, 'tsv', 'tsv')
        ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

        // Clean up scrubbed reads
        QC(SCRUBBER.out.scrubbed_extra)
        ch_versions = ch_versions.mix(QC.out.versions.first())
    } else {
        // Clean up raw reads
        QC(GATHER.out.raw_fastq)
        ch_versions = ch_versions.mix(QC.out.versions.first())
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
