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
if (params.wf == 'kleborate') include { KLEBORATE } from '../subworkflows/local/kleborate/main';
if (params.wf == 'mashtree') include { SPATYPER } from '../subworkflows/local/mashtree/main';
if (params.wf == 'spatyper') include { SPATYPER } from '../subworkflows/local/spatyper/main';
if (params.wf == 'staphtyper') include { STAPHTYPER } from '../subworkflows/local/staphtyper/main';
if (params.wf == 'staphopiasccmec') include { STAPHOPIASCCMEC } from '../subworkflows/local/staphopiasccmec/main';

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' addParams( options: [publish_to_base: true] )

/*
include { FASTANI } from '../modules/nf-core/modules/fastani/main'
include { GTDBTK_CLASSIFYWF } from '../modules/nf-core/modules/gtdbtk/classifywf/main'
include { HICAP } from '../modules/nf-core/modules/hicap/main'
include { IQTREE } from '../modules/nf-core/modules/iqtree/main'
include { ISMAPPER } from '../modules/nf-core/modules/ismapper/main'
include { PIRATE } from '../modules/nf-core/modules/pirate/main'
include { PROKKA } from '../modules/nf-core/modules/prokka/main'
include { ROARY } from '../modules/nf-core/modules/roary/main'
include { SNPDISTS } from '../modules/nf-core/modules/snpdists/main'
include { SPATYPER } from '../modules/nf-core/modules/spatyper/main'
*/
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
    } else if (params.wf == 'kleborate') {
        KLEBORATE(samples)
        ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    }  else if (params.wf == 'spatyper') {
        SPATYPER(samples)
        ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    }  else if (params.wf == 'staphtyper') {
        STAPHTYPER(samples)
        ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    } else if (params.wf == 'staphopiasccmec') {
        STAPHOPIASCCMEC(samples)
        ch_versions = ch_versions.mix(AGRVATE.out.versions)
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
