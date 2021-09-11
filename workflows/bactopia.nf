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
include { create_input_channel; setup_datasets } from '../modules/local/bactopia/inputs.nf'
include { get_resources; print_efficiency } from '../modules/local/utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/


// Core
include { ANNOTATE_GENOME } from '../modules/local/bactopia/annotate_genome/main.nf'
include { ASSEMBLE_GENOME } from '../modules/local/bactopia/assemble_genome/main.nf'
include { ASSEMBLY_QC } from '../modules/local/bactopia/assembly_qc/main.nf'
include { GATHER_SAMPLES } from '../modules/local/bactopia/gather_samples/main.nf'
include { MINMER_SKETCH } from '../modules/local/bactopia/minmer_sketch/main.nf'
include { QC_READS } from '../modules/local/bactopia/qc_reads/main.nf'

// Require Datasets
include { ANTIMICROBIAL_RESISTANCE } from '../modules/local/bactopia/antimicrobial_resistance/main.nf'
include { ARIBA_ANALYSIS } from '../modules/local/bactopia/ariba_analysis/main.nf'
include { BLAST } from '../modules/local/bactopia/blast/main.nf'
include { CALL_VARIANTS } from '../modules/local/bactopia/call_variants/main.nf'
include { MAPPING_QUERY } from '../modules/local/bactopia/mapping_query/main.nf'
include { MINMER_QUERY } from '../modules/local/bactopia/minmer_query/main.nf'
include { SEQUENCE_TYPE } from '../modules/local/bactopia/sequence_type/main.nf'

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
    print_efficiency(RESOURCES.MAX_CPUS) 
    datasets = setup_datasets()

    // Core steps
    GATHER_SAMPLES(create_input_channel(run_type, datasets['genome_size']))
    QC_READS(GATHER_SAMPLES.out.raw_fastq)
    ASSEMBLE_GENOME(QC_READS.out.fastq_assembly)
    ASSEMBLY_QC(ASSEMBLE_GENOME.out.fna, Channel.fromList(['checkm', 'quast']))
    ANNOTATE_GENOME(ASSEMBLE_GENOME.out.fna, Channel.fromPath(datasets['proteins']), Channel.fromPath(datasets['training_set']))
    MINMER_SKETCH(QC_READS.out.fastq)

    // Optional steps that require datasets
    // Species agnostic
    ANTIMICROBIAL_RESISTANCE(ANNOTATE_GENOME.out.annotations, datasets['amr'])
    ARIBA_ANALYSIS(QC_READS.out.fastq, datasets['ariba'])
    MINMER_QUERY(MINMER_SKETCH.out.sketch, datasets['minmer'])

    // Species Specific
    BLAST(ASSEMBLE_GENOME.out.blastdb, datasets['blast'])
    CALL_VARIANTS(QC_READS.out.fastq, datasets['references'])
    MAPPING_QUERY(QC_READS.out.fastq, datasets['mapping'])
    SEQUENCE_TYPE(ASSEMBLE_GENOME.out.fna_fastq, datasets['mlst'])
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
