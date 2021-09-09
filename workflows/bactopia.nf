#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
WorkflowMain.initialise(workflow, params, log, schema_filename='conf/schema/bactopia.json')
run_type = WorkflowBactopia.initialise(workflow, params, log)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
include { GATHER_SAMPLES } from '../modules/local/bactopia/gather_samples/main.nf'
include { QC_READS } from '../modules/local/bactopia/qc_reads/main.nf'
include { ASSEMBLE_GENOME } from '../modules/local/bactopia/assemble_genome/main.nf'
include { ASSEMBLY_QC } from '../modules/local/bactopia/assembly_qc/main.nf'
include { ANNOTATE_GENOME } from '../modules/local/bactopia/annotate_genome/main.nf'
include { SEQUENCE_TYPE } from '../modules/local/bactopia/sequence_type/main.nf'
include { ARIBA_ANALYSIS } from '../modules/local/bactopia/ariba_analysis/main.nf'
include { MINMER_SKETCH } from '../modules/local/bactopia/minmer_sketch/main.nf'
include { MINMER_QUERY } from '../modules/local/bactopia/minmer_query/main.nf'
include { CALL_VARIANTS } from '../modules/local/bactopia/call_variants/main.nf'
include { ANTIMICROBIAL_RESISTANCE } from '../modules/local/bactopia/antimicrobial_resistance/main.nf'
include { BLAST } from '../modules/local/bactopia/blast/main.nf'
include { MAPPING_QUERY } from '../modules/local/bactopia/mapping_query/main.nf'
include { create_input_channel } from '../modules/local/utilities/bactopia-channels.nf'
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

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow BACTOPIA {
    SPECIES_GENOME_SIZE = params.genome_size
    PROKKA_PROTEINS = params.empty_proteins
    PRODIGAL_TF = params.empty_tf
    GATHER_SAMPLES(create_input_channel(run_type, params.genome_size))
    QC_READS(GATHER_SAMPLES.out.raw_fastq)
    ASSEMBLE_GENOME(QC_READS.out.fastq_assembly)
    ANNOTATE_GENOME(ASSEMBLE_GENOME.out.fna, Channel.fromPath(params.empty_proteins), Channel.fromPath(params.empty_tf))
    MINMER_SKETCH(QC_READS.out.fastq)
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
    Command Line    : ${workflow.commandLine}
    Resumed         : ${workflow.resume}
    Completed At    : ${workflow.complete}
    Duration        : ${workflow.duration}
    Success         : ${workflow.success}
    Exit Code       : ${workflow.exitStatus}
    Error Report    : ${workflow.errorReport ?: '-'}
    Launch Dir      : ${workflow.launchDir}
    """
}

/*
========================================================================================
    THE END
========================================================================================
*/
