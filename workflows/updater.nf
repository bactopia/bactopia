#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    NOTE TO USER

    This Nextflow script is really only meant to be used by Bactopia developers. Now we
    can't stop you from using it, but before you do consider if its needed. The datbases
    built in the workflow are already available to you in Bactopia. If there is truely
    a need to run this workflow, please consider opening an issue on GitHub so we can
    discuss it. If you see a need, you likely are not alone and maybe it warrants a
    version update for Bactopia!
========================================================================================
*/

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
include { get_resources; get_schemas; print_efficiency } from '../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { AMRFINDERPLUS_UPDATE } from '../modules/nf-core/amrfinderplus/update/main'
include { MLST_UPDATE } from '../modules/nf-core/mlst/update/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'  addParams( options: [publish_to_base: true] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow UPDATER {
    ch_versions = Channel.empty()

    AMRFINDERPLUS_UPDATE()
    ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions.first())

    MLST_UPDATE()
    ch_versions = ch_versions.mix(MLST_UPDATE.out.versions.first())

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
    Bactopia Updater Execution Summary
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
