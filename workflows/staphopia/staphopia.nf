#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
include { create_input_channel; check_input_fofn; } from '../lib/nf/bactopia'
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

include { AMRFINDERPLUS } from '../subworkflows/local/amrfinderplus/main'
include { ASSEMBLER } from '../subworkflows/local/assembler/main'
include { DATASETS } from '../modules/local/bactopia/datasets/main'
include { GATHER } from '../subworkflows/local/gather/main'
include { SKETCHER } from '../subworkflows/local/sketcher/main'
include { MLST } from '../subworkflows/local/mlst/main'
include { QC } from '../subworkflows/local/qc/main'


// Annotation wih Bakta or Prokka
if (params.use_bakta) {
    include { BAKTA as ANNOTATOR } from '../subworkflows/local/bakta/main'
} else {
    include { PROKKA as ANNOTATOR } from '../subworkflows/local/prokka/main'
}

// Staphylococcus aureus specific
include { STAPHTYPER } from '../subworkflows/local/staphtyper/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS as DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main' addParams( options: [process_name:"bactopia"] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow STAPHOPIA {
    ch_versions = Channel.empty()

    // Core steps
    DATASETS()
    GATHER(create_input_channel(run_type, params.genome_size, params.species))
    QC(GATHER.out.raw_fastq)
    ASSEMBLER(QC.out.fastq)
    SKETCHER(ASSEMBLER.out.fna, DATASETS.out.mash_db, DATASETS.out.sourmash_db)
    ANNOTATOR(ASSEMBLER.out.fna)
    AMRFINDERPLUS(ANNOTATOR.out.annotations, DATASETS.out.amrfinderplus_db)
    MLST(ASSEMBLER.out.fna, DATASETS.out.mlst_db)

    // Staphylococcus aureus specific
    STAPHTYPER(ASSEMBLER.out.fna)

    // Collect Versions
    ch_versions = ch_versions.mix(GATHER.out.versions.first())
    ch_versions = ch_versions.mix(QC.out.versions.first())
    ch_versions = ch_versions.mix(ASSEMBLER.out.versions.first())
    ch_versions = ch_versions.mix(ANNOTATOR.out.versions.first())
    ch_versions = ch_versions.mix(SKETCHER.out.versions.first())
    ch_versions = ch_versions.mix(AMRFINDERPLUS.out.versions.first())
    ch_versions = ch_versions.mix(MLST.out.versions.first())
    ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile())
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
    Staphopia Execution Summary
    -------------------------------
    Workflow         : ${params.wf}
    Bactopia Version : ${workflow.manifest.version}
    Nextflow Version : ${nextflow.version}
    Command Line     : ${workflow.commandLine}
    Profile          : ${workflow.profile}
    Completed At     : ${workflow.complete}
    Duration         : ${workflow.duration}
    Success          : ${workflow.success}
    Exit Code        : ${workflow.exitStatus}
    Error Report     : ${workflow.errorReport ?: '-'}
    Launch Dir       : ${workflow.launchDir}
    ${colors.bgreen}Merged Results${colors.reset}   : ${colors.green}${params.outdir}/bactopia-runs/${params.rundir}${colors.reset}

    Further analyze your samples using Bactopia Tools, with the following command:
    --------------------------------------------------------------------------------
    ${colors.cyan}bactopia -profile ${workflow.profile} --bactopia ${params.outdir} --wf <REPLACE_WITH_BACTOPIA_TOOL_NAME>${colors.reset}

    Examples:
    ${colors.cyan}bactopia -profile ${workflow.profile} --bactopia ${params.outdir} --wf pangenome${colors.reset}
    ${colors.cyan}bactopia -profile ${workflow.profile} --bactopia ${params.outdir} --wf mobsuite${colors.reset}
    ${colors.cyan}bactopia -profile ${workflow.profile} --bactopia ${params.outdir} --wf busco${colors.reset}

    See the full list of available Bactopia Tools: ${colors.cyan}bactopia --list_wfs${colors.reset}
    """
}

/*
========================================================================================
    THE END
========================================================================================
*/
