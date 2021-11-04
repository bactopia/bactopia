#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
include { collect_samples } from '../lib/nf/bactopia_tools'
include { get_resources; get_schemas; print_efficiency } from '../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
SCHEMAS = get_schemas()
WorkflowMain.initialise(workflow, params, log, schema_filename=SCHEMAS)
WorkflowBactopiaTools.initialise(workflow, params, log, schema_filename=SCHEMAS)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/* All certain steps to be rerun
include { ANNOTATE_GENOME } from '../modules/local/bactopia/annotate_genome/main'
include { ASSEMBLE_GENOME } from '../modules/local/bactopia/assemble_genome/main'
include { ANTIMICROBIAL_RESISTANCE } from '../modules/local/bactopia/antimicrobial_resistance/main'
include { BLAST } from '../modules/local/bactopia/blast/main'
include { CALL_VARIANTS } from '../modules/local/bactopia/call_variants/main'
include { MAPPING_QUERY } from '../modules/local/bactopia/mapping_query/main'
include { MINMER_QUERY } from '../modules/local/bactopia/minmer_query/main'
*/

// Subworkflows
if (params.wf == 'agrvate') include { AGRVATE } from '../subworkflows/local/agrvate/main';
if (params.wf == 'ectyper') include { ECTYPER } from '../subworkflows/local/ectyper/main';
//if (params.wf == 'eggnog') include { EGGNOG } from '../subworkflows/local/eggnog/main';
//if (params.wf == 'emmtyper') include { EMMTYPER } from '../subworkflows/local/emmtyper/main';
//if (params.wf == 'fastani') include { FASTANI } from '../subworkflows/local/fastani/main';
//if (params.wf == 'gtdb') include { GTDB } from '../subworkflows/local/gtdb/main';
if (params.wf == 'hicap') include { HICAP } from '../subworkflows/local/hicap/main';
if (params.wf == 'kleborate') include { KLEBORATE } from '../subworkflows/local/kleborate/main';
if (params.wf == 'mashtree') include { MASHTREE } from '../subworkflows/local/mashtree/main';
//if (params.wf == 'meningotype') include { MENINGOTYPE } from '../subworkflows/local/meningotype/main';
//if (params.wf == 'ngmaster') include { NGMASTER } from '../subworkflows/local/ngmaster/main';
//if (params.wf == 'seqsero2') include { SEQSERO2 } from '../subworkflows/local/seqsero2/main';
if (params.wf == 'spatyper') include { SPATYPER } from '../subworkflows/local/spatyper/main';
if (params.wf == 'staphtyper') include { STAPHTYPER } from '../subworkflows/local/staphtyper/main';
if (params.wf == 'staphopiasccmec') include { STAPHOPIASCCMEC } from '../subworkflows/local/staphopiasccmec/main';
//if (params.wf == 'tbprofiler') include { SEQSERO2 } from '../subworkflows/local/tbprofiler/main';

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' addParams( options: [publish_to_base: true] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow BACTOPIATOOLS {
    print_efficiency(RESOURCES.MAX_CPUS)
    samples = Channel.fromList(collect_samples(params.bactopia, params.workflows[params.wf].ext, params.include, params.exclude))
    ch_versions = Channel.empty()

    if (params.wf == 'agrvate') {
        AGRVATE(samples)
        ch_versions = ch_versions.mix(AGRVATE.out.versions)
    } else if (params.wf == 'hicap') {
        HICAP(samples)
        ch_versions = ch_versions.mix(HICAP.out.versions)
    }  else if (params.wf == 'kleborate') {
        KLEBORATE(samples)
        ch_versions = ch_versions.mix(KLEBORATE.out.versions)
    } else if (params.wf == 'mashtree') {
        MASHTREE(samples)
        ch_versions = ch_versions.mix(MASHTREE.out.versions)
    } else if (params.wf == 'spatyper') {
        SPATYPER(samples)
        ch_versions = ch_versions.mix(SPATYPER.out.versions)
    } else if (params.wf == 'staphtyper') {
        STAPHTYPER(samples)
        ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    } else if (params.wf == 'staphopiasccmec') {
        STAPHOPIASCCMEC(samples)
        ch_versions = ch_versions.mix(STAPHOPIASCCMEC.out.versions)
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
    Bactopia Tools: `${params.wf} Execution Summary
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
