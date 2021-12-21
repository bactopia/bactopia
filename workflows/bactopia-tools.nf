#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
include { collect_samples; collect_local_files } from '../lib/nf/bactopia_tools'
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
if (params.containsKey('accession')) {
    // Only import when made available
    include { NCBIGENOMEDOWNLOAD } from '../subworkflows/local/ncbigenomedownload/main'
}
if (params.wf == 'agrvate') include { AGRVATE } from '../subworkflows/local/agrvate/main';
if (params.wf == 'bakta') include { BAKTA } from '../subworkflows/local/bakta/main';
if (params.wf == 'ectyper') include { ECTYPER } from '../subworkflows/local/ectyper/main';
if (params.wf == 'eggnog') include { EGGNOG } from '../subworkflows/local/eggnog/main';
if (params.wf == 'emmtyper') include { EMMTYPER } from '../subworkflows/local/emmtyper/main';
if (params.wf == 'fastani') include { FASTANI } from '../subworkflows/local/fastani/main';
if (params.wf == 'gtdb') include { GTDB } from '../subworkflows/local/gtdb/main';
if (params.wf == 'hicap') include { HICAP } from '../subworkflows/local/hicap/main';
if (params.wf == 'kleborate') include { KLEBORATE } from '../subworkflows/local/kleborate/main';
if (params.wf == 'lissero') include { LISSERO } from '../subworkflows/local/lissero/main';
if (params.wf == 'mashtree') include { MASHTREE } from '../subworkflows/local/mashtree/main';
if (params.wf == 'meningotype') include { MENINGOTYPE } from '../subworkflows/local/meningotype/main';
if (params.wf == 'ngmaster') include { NGMASTER } from '../subworkflows/local/ngmaster/main';
if (params.wf == 'pangenome') include { PANGENOME } from '../subworkflows/local/pangenome/main';
if (params.wf == 'pangenome') include { PROKKA } from '../subworkflows/local/prokka/main';
if (params.wf == 'seqsero2') include { SEQSERO2 } from '../subworkflows/local/seqsero2/main';
if (params.wf == 'spatyper') include { SPATYPER } from '../subworkflows/local/spatyper/main';
if (params.wf == 'staphtyper') include { STAPHTYPER } from '../subworkflows/local/staphtyper/main';
if (params.wf == 'staphopiasccmec') include { STAPHOPIASCCMEC } from '../subworkflows/local/staphopiasccmec/main';
if (params.wf == 'tbprofiler') include { TBPROFILER } from '../subworkflows/local/tbprofiler/main';

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
    ch_versions = Channel.empty()
    ch_local_samples = Channel.fromList(collect_samples(params.bactopia, params.workflows[params.wf].ext, params.include, params.exclude))
    ch_local_files = Channel.fromList(collect_local_files(params.containsKey('assembly') ? params.assembly : null , params.containsKey('assembly_pattern') ? params.assembly_pattern : null))
    
    // Include public genomes (optional)
    ch_gather_files = Channel.empty()
    if (params.containsKey('accession')) {
        ch_downloads = Channel.empty()
        if (params.accession || params.accessions || params.species) {
            NCBIGENOMEDOWNLOAD()
            ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD.out.versions.first())
            ch_downloads = NCBIGENOMEDOWNLOAD.out.bactopia_tools.mix(ch_local_files)
        } else {
            ch_downloads = ch_local_files
        }

        if (params.wf == 'pangenome') {
            // Create prokka gff3 files for PIRATE/Roary
            PROKKA(ch_downloads)
            ch_versions = ch_versions.mix(PROKKA.out.versions.first())
            ch_gather_files = PROKKA.out.gff
        } else {
            ch_gather_files = ch_downloads
        }
    } else {
        ch_gather_files = ch_local_files
    }

    // Include local GFF files (optional)
    ch_final_downloads = ch_gather_files.mix(Channel.fromList(collect_local_files(params.containsKey('gff') ? params.gff : null , params.containsKey('gff_pattern') ? params.gff_pattern : null)))

    samples = ch_local_samples.mix(ch_final_downloads)
    if (params.wf == 'agrvate') {
        AGRVATE(samples)
        ch_versions = ch_versions.mix(AGRVATE.out.versions)
    } else if (params.wf == 'bakta') {
        BAKTA(samples)
        ch_versions = ch_versions.mix(BAKTA.out.versions)
    } else if (params.wf == 'ectyper') {
        ECTYPER(samples)
        ch_versions = ch_versions.mix(ECTYPER.out.versions)
    } else if (params.wf == 'eggnog') {
        EGGNOG(samples)
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
    } else if (params.wf == 'emmtyper') {
        EMMTYPER(samples)
        ch_versions = ch_versions.mix(EMMTYPER.out.versions)
    } else if (params.wf == 'fastani') {
        FASTANI(samples, ch_downloads)
        ch_versions = ch_versions.mix(FASTANI.out.versions)
    } else if (params.wf == 'gtdb') {
        GTDB(samples)
        ch_versions = ch_versions.mix(GTDB.out.versions)
    }  else if (params.wf == 'hicap') {
        HICAP(samples)
        ch_versions = ch_versions.mix(HICAP.out.versions)
    } else if (params.wf == 'kleborate') {
        KLEBORATE(samples)
        ch_versions = ch_versions.mix(KLEBORATE.out.versions)
    } else if (params.wf == 'lissero') {
        LISSERO(samples)
        ch_versions = ch_versions.mix(LISSERO.out.versions)
    } else if (params.wf == 'mashtree') {
        MASHTREE(samples)
        ch_versions = ch_versions.mix(MASHTREE.out.versions)
    } else if (params.wf == 'meningotype') {
        MENINGOTYPE(samples)
        ch_versions = ch_versions.mix(MENINGOTYPE.out.versions)
    } else if (params.wf == 'ngmaster') {
        NGMASTER(samples)
        ch_versions = ch_versions.mix(NGMASTER.out.versions)
    } else if (params.wf == 'pangenome') {
        samples.collect{meta, gff -> gff}.map{ gff -> [[id: params.use_roary ? 'roary' : 'pirate'], gff]}.set{ ch_merge_gff }
        PANGENOME(ch_merge_gff)
        ch_versions = ch_versions.mix(PANGENOME.out.versions)
    } else if (params.wf == 'seqsero2') {
        SEQSERO2(samples)
        ch_versions = ch_versions.mix(SEQSERO2.out.versions)
    } else if (params.wf == 'spatyper') {
        SPATYPER(samples)
        ch_versions = ch_versions.mix(SPATYPER.out.versions)
    } else if (params.wf == 'staphtyper') {
        STAPHTYPER(samples)
        ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    } else if (params.wf == 'staphopiasccmec') {
        STAPHOPIASCCMEC(samples)
        ch_versions = ch_versions.mix(STAPHOPIASCCMEC.out.versions)
    } else if (params.wf == 'tbprofiler') {
        TBPROFILER(samples)
        ch_versions = ch_versions.mix(TBPROFILER.out.versions)
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
