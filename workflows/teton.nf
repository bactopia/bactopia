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
include { BACTOPIA_SAMPLESHEET } from '../modules/local/bactopia/teton/main'
include { BRACKEN } from '../subworkflows/local/bracken/main'
include { SCRUBBER } from '../subworkflows/local/scrubber/main' addParams( options: [ignore: [".fna.gz"]] )
include { CSVTK_JOIN } from '../modules/nf-core/csvtk/join/main' addParams( options: [process_name: "teton-report"] )
include { CSVTK_CONCAT } from '../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'teton-concat', process_name: params.merge_folder] )
include { CSVTK_CONCAT as CSVTK_CONCAT_SIZEMEUP } from '../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'sizemeup', process_name: params.merge_folder] )
include { CSVTK_CONCAT as CSVTK_CONCAT_BACTERIA } from '../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'teton-prepare-bacteria', process_name: params.merge_folder] )
include { CSVTK_CONCAT as CSVTK_CONCAT_NONBACTERIA } from '../modules/nf-core/csvtk/concat/main' addParams( options: [logs_subdir: 'teton-prepare-nonbacteria', process_name: params.merge_folder] )

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
    ch_versions = Channel.empty()
    ch_merged_teton = Channel.empty()
    ch_join_teton = Channel.empty()

    // Core Steps
    GATHER(create_input_channel(run_type, params.genome_size, params.species))
    ch_versions = ch_versions.mix(GATHER.out.versions.first())

    // Remove host reads
    SCRUBBER(GATHER.out.fastq_only)
    ch_versions = ch_versions.mix(SCRUBBER.out.versions)

    // Taxon Classification & Abundance
    BRACKEN(SCRUBBER.out.scrubbed)
    ch_versions = ch_versions.mix(BRACKEN.out.versions)

    // Determine genome size and create sample sheet
    BACTOPIA_SAMPLESHEET(BRACKEN.out.teton_classification)
    ch_versions = ch_versions.mix(BACTOPIA_SAMPLESHEET.out.versions)

    // Join Scrubber and Bracken results
    CSVTK_JOIN(SCRUBBER.out.tsv.join(BRACKEN.out.tsv, by:[0]), 'tsv', 'tsv', 'sample')
    CSVTK_JOIN.out.csv.collect{meta, csv -> csv}.map{ csv -> [[id:'teton'], csv]}.set{ ch_merge_teton }
    ch_versions = ch_versions.mix(CSVTK_JOIN.out.versions)

    // Merge the results
    CSVTK_CONCAT(ch_merge_teton, 'tsv', 'tsv')
    ch_merged_teton = ch_merged_teton.mix(CSVTK_CONCAT.out.csv)
    BACTOPIA_SAMPLESHEET.out.bacteria_tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'teton-prepare'], tsv]}.set{ ch_merge_prepare }
    CSVTK_CONCAT_BACTERIA(ch_merge_prepare, 'tsv', 'tsv')
    BACTOPIA_SAMPLESHEET.out.nonbacteria_tsv.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'teton-prepare-nonbacteria'], tsv]}.set{ ch_merge_prepare_non }
    CSVTK_CONCAT_NONBACTERIA(ch_merge_prepare_non, 'tsv', 'tsv')
    BACTOPIA_SAMPLESHEET.out.sizemeup.collect{meta, tsv -> tsv}.map{ tsv -> [[id:'sizemeup'], tsv]}.set{ ch_merge_sizemeup }
    CSVTK_CONCAT_SIZEMEUP(ch_merge_sizemeup, 'tsv', 'tsv')
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
    def colors = NfcoreTemplate.logColours(params.monochrome_logs)

    println """
    Teton Execution Summary
    ---------------------------
    Workflow         : ${params.wf}
    Bactopia Version : ${workflow.manifest.version}
    Nextflow Version : ${nextflow.version}
    Command Line     : ${workflow.commandLine}
    Profile          : ${workflow.profile}
    Resumed          : ${workflow.resume}
    Completed At     : ${workflow.complete}
    Duration         : ${workflow.duration}
    Success          : ${workflow.success}
    Exit Code        : ${workflow.exitStatus}
    Error Report     : ${workflow.errorReport ?: '-'}
    Launch Dir       : ${workflow.launchDir}
    ${colors.bgreen}Merged Results${colors.reset}   : ${colors.green}${params.outdir}/bactopia-runs/${params.rundir}${colors.reset}
    
    Further analyze bacterial samples using Bactopia, with the following command:
    --------------------------------------------------------------------------------
    ${colors.cyan}bactopia -profile ${workflow.profile} --samples ${params.outdir}/bactopia-runs/${params.rundir}/merged-results/teton-prepare.tsv${colors.reset}
    """
}

/*
========================================================================================
    LA FIN
========================================================================================
*/
