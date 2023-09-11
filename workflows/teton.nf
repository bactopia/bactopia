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
include { BRACKEN } from '../subworkflows/local/bracken/main'
include { MIDAS } from '../subworkflows/local/midas/main'
include { SCRUBBER } from '../subworkflows/local/scrubber/main' addParams( options: [publish_to_base: true, ignore: [".fna.gz"]] )
include { CSVTK_JOIN } from '../modules/nf-core/csvtk/join/main' addParams( options: [publish_to_base: true] )
include { CSVTK_CONCAT } from '../modules/nf-core/csvtk/concat/main' addParams( options: [publish_to_base: true] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main' addParams( options: [publish_to_base: true] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow TETON {
    print_efficiency(RESOURCES.MAX_CPUS)
    ch_versions = Channel.empty()
    ch_merged_teton = Channel.empty()
    ch_join_teton = Channel.empty()

    // Core Steps
    GATHER_SAMPLES(create_input_channel(run_type, params.genome_size))
    ch_versions = ch_versions.mix(GATHER_SAMPLES.out.versions.first())

    // Remove host reads
    SCRUBBER(GATHER_SAMPLES.out.raw_fastq)
    ch_versions = ch_versions.mix(SCRUBBER.out.versions)

    // Taxon Classification & Abundance
    BRACKEN(SCRUBBER.out.scrubbed)
    ch_versions = ch_versions.mix(BRACKEN.out.versions)

    if (params.midas_db) {
        // Species Abundance
        MIDAS(SCRUBBER.out.scrubbed)
        ch_versions = ch_versions.mix(MIDAS.out.versions)

        // Join Bracken and MIDAS results
        CSVTK_JOIN(BRACKEN.out.tsv.join(MIDAS.out.tsv, by:[0]), 'tsv', 'tsv', 'sample')
        CSVTK_JOIN.out.csv.collect{meta, csv -> csv}.map{ csv -> [[id:'teton'], csv]}.set{ ch_merge_teton }
        ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)
    } else {
        BRACKEN.out.tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'teton'], tsv]}.set{ ch_merge_teton }
    }

    // Join the results
    CSVTK_CONCAT(ch_merge_teton, 'tsv', 'tsv')
    ch_merged_teton = ch_merged_teton.mix(CSVTK_CONCAT.out.csv)
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)

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
    LA FIN
========================================================================================
*/
