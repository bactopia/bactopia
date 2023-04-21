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
include { BLAST } from '../modules/local/bactopia/blast/main'
include { MAPPING_QUERY } from '../modules/local/bactopia/mapping_query/main'
include { MINMER_QUERY } from '../modules/local/bactopia/minmer_query/main'
*/

// Subworkflows
if (params.containsKey('accession')) {
    // Only import when made available
    include { NCBIGENOMEDOWNLOAD } from '../subworkflows/local/ncbigenomedownload/main'
}
if (params.wf == 'abricate') include { ABRICATE } from '../subworkflows/local/abricate/main';
if (params.wf == 'agrvate') include { AGRVATE } from '../subworkflows/local/agrvate/main';
if (params.wf == 'amrfinderplus') include { AMRFINDERPLUS } from '../subworkflows/local/amrfinderplus/main';
if (params.wf == 'ariba') include { ARIBA } from '../subworkflows/local/ariba/main';
if (params.wf == 'bakta') include { BAKTA } from '../subworkflows/local/bakta/main';
if (params.wf == 'blastn') include { BLASTN } from '../subworkflows/local/blastn/main';
if (params.wf == 'blastp') include { BLASTP } from '../subworkflows/local/blastp/main';
if (params.wf == 'blastx') include { BLASTX } from '../subworkflows/local/blastx/main';
if (params.wf == 'bracken') include { BRACKEN } from '../subworkflows/local/bracken/main';
if (params.wf == 'busco') include { BUSCO } from '../subworkflows/local/busco/main';
if (params.wf == 'checkm') include { CHECKM } from '../subworkflows/local/checkm/main';
if (params.wf == 'ectyper') include { ECTYPER } from '../subworkflows/local/ectyper/main';
if (params.wf == 'eggnog') include { EGGNOG } from '../subworkflows/local/eggnog/main';
if (params.wf == 'emmtyper') include { EMMTYPER } from '../subworkflows/local/emmtyper/main';
if (params.wf == 'fastani') include { FASTANI } from '../subworkflows/local/fastani/main';
if (params.wf == 'gamma') include { GAMMA } from '../subworkflows/local/gamma/main';
if (params.wf == 'genotyphi') include { GENOTYPHI } from '../subworkflows/local/genotyphi/main';
if (params.wf == 'gtdb') include { GTDB } from '../subworkflows/local/gtdb/main';
if (params.wf == 'hicap') include { HICAP } from '../subworkflows/local/hicap/main';
if (params.wf == 'hpsuissero') include { HPSUISSERO } from '../subworkflows/local/hpsuissero/main';
if (params.wf == 'ismapper') include { ISMAPPER } from '../subworkflows/local/ismapper/main';
if (params.wf == 'kleborate') include { KLEBORATE } from '../subworkflows/local/kleborate/main';
if (params.wf == 'kraken2') include { KRAKEN2 } from '../subworkflows/local/kraken2/main';
if (params.wf == 'legsta') include { LEGSTA } from '../subworkflows/local/legsta/main';
if (params.wf == 'lissero') include { LISSERO } from '../subworkflows/local/lissero/main';
if (params.wf == 'mashdist') include { MASHDIST } from '../subworkflows/local/mashdist/main';
if (params.wf == 'mashtree') include { MASHTREE } from '../subworkflows/local/mashtree/main';
if (params.wf == 'mcroni') include { MCRONI } from '../subworkflows/local/mcroni/main';
if (params.wf == 'meningotype') include { MENINGOTYPE } from '../subworkflows/local/meningotype/main';
if (params.wf == 'merlin') include { MERLIN } from '../subworkflows/local/merlin/main';
if (params.wf == 'midas') include { MIDAS } from '../subworkflows/local/midas/main';
if (params.wf == 'mlst') include { MLST } from '../subworkflows/local/mlst/main';
if (params.wf == 'mobsuite') include { MOBSUITE } from '../subworkflows/local/mobsuite/main';
if (params.wf == 'mykrobe') include { MYKROBE } from '../subworkflows/local/mykrobe/main';
if (params.wf == 'ngmaster') include { NGMASTER } from '../subworkflows/local/ngmaster/main';
if (params.wf == 'pangenome') include { PANGENOME } from '../subworkflows/local/pangenome/main';
if (params.wf == 'pangenome') include { PROKKA } from '../subworkflows/local/prokka/main';
if (params.wf == 'pasty') include { PASTY } from '../subworkflows/local/pasty/main';
if (params.wf == 'pbptyper') include { PBPTYPER } from '../subworkflows/local/pbptyper/main';
if (params.wf == 'plasmidfinder') include { PLASMIDFINDER } from '../subworkflows/local/plasmidfinder/main';
if (params.wf == 'rgi') include { RGI } from '../subworkflows/local/rgi/main';
if (params.wf == 'seqsero2') include { SEQSERO2 } from '../subworkflows/local/seqsero2/main';
if (params.wf == 'seroba') include { SEROBA } from '../subworkflows/local/seroba/main';
if (params.wf == 'shigatyper') include { SHIGATYPER } from '../subworkflows/local/shigatyper/main';
if (params.wf == 'shigeifinder') include { SHIGEIFINDER } from '../subworkflows/local/shigeifinder/main';
if (params.wf == 'sistr') include { SISTR } from '../subworkflows/local/sistr/main';
if (params.wf == 'snippy') include { SNIPPY } from '../subworkflows/local/snippy/main';
if (params.wf == 'spatyper') include { SPATYPER } from '../subworkflows/local/spatyper/main';
if (params.wf == 'ssuissero') include { SSUISSERO } from '../subworkflows/local/ssuissero/main';
if (params.wf == 'staphtyper') include { STAPHTYPER } from '../subworkflows/local/staphtyper/main';
if (params.wf == 'staphopiasccmec') include { STAPHOPIASCCMEC } from '../subworkflows/local/staphopiasccmec/main';
if (params.wf == 'stecfinder') include { STECFINDER } from '../subworkflows/local/stecfinder/main';
if (params.wf == 'tbprofiler') include { TBPROFILER } from '../subworkflows/local/tbprofiler/main';
if (params.wf == 'tblastn') include { TBLASTN } from '../subworkflows/local/tblastn/main';
if (params.wf == 'tblastx') include { TBLASTX } from '../subworkflows/local/tblastx/main';
if (params.wf == 'teton') include { TETON } from '../subworkflows/local/teton/main';

if (['amrfinderplus', 'mlst'].contains(params.wf)) {
    include { DATASETS } from '../modules/local/bactopia/datasets/main'
}

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

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
    if (params.wf == 'abricate') {
        ABRICATE(samples)
        ch_versions = ch_versions.mix(ABRICATE.out.versions)
    } else if (params.wf == 'agrvate') {
        AGRVATE(samples)
        ch_versions = ch_versions.mix(AGRVATE.out.versions)
    } else if (params.wf == 'amrfinderplus') {
        DATASETS()
        AMRFINDERPLUS(samples, DATASETS.out.amrfinderplus_db)
        ch_versions = ch_versions.mix(AMRFINDERPLUS.out.versions)
    } else if (params.wf == 'ariba') {
        ARIBA(samples)
        ch_versions = ch_versions.mix(ARIBA.out.versions)
    } else if (params.wf == 'bakta') {
        BAKTA(samples)
        ch_versions = ch_versions.mix(BAKTA.out.versions)
    } else if (params.wf == 'blastn') {
        BLASTN(samples)
        ch_versions = ch_versions.mix(BLASTN.out.versions)
    } else if (params.wf == 'blastp') {
        BLASTP(samples)
        ch_versions = ch_versions.mix(BLASTP.out.versions)
    } else if (params.wf == 'blastx') {
        BLASTX(samples)
        ch_versions = ch_versions.mix(BLASTX.out.versions)
    } else if (params.wf == 'bracken') {
        BRACKEN(samples)
        ch_versions = ch_versions.mix(BRACKEN.out.versions)
    } else if (params.wf == 'busco') {
        BUSCO(samples)
        ch_versions = ch_versions.mix(BUSCO.out.versions)
    } else if (params.wf == 'checkm') {
        CHECKM(samples)
        ch_versions = ch_versions.mix(CHECKM.out.versions)
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
    } else if (params.wf == 'gamma') {
        GAMMA(samples)
        ch_versions = ch_versions.mix(GAMMA.out.versions)
    } else if (params.wf == 'genotyphi') {
        GENOTYPHI(samples)
        ch_versions = ch_versions.mix(GENOTYPHI.out.versions)
    } else if (params.wf == 'gtdb') {
        GTDB(samples)
        ch_versions = ch_versions.mix(GTDB.out.versions)
    } else if (params.wf == 'hicap') {
        HICAP(samples)
        ch_versions = ch_versions.mix(HICAP.out.versions)
    } else if (params.wf == 'hpsuissero') {
        HPSUISSERO(samples)
        ch_versions = ch_versions.mix(HPSUISSERO.out.versions)
    } else if (params.wf == 'ismapper') {
        ISMAPPER(samples)
        ch_versions = ch_versions.mix(ISMAPPER.out.versions)
    } else if (params.wf == 'kleborate') {
        KLEBORATE(samples)
        ch_versions = ch_versions.mix(KLEBORATE.out.versions)
    } else if (params.wf == 'kraken2') {
        KRAKEN2(samples)
        ch_versions = ch_versions.mix(KRAKEN2.out.versions)
    } else if (params.wf == 'legsta') {
        LEGSTA(samples)
        ch_versions = ch_versions.mix(LEGSTA.out.versions)
    } else if (params.wf == 'lissero') {
        LISSERO(samples)
        ch_versions = ch_versions.mix(LISSERO.out.versions)
    } else if (params.wf == 'mashdist') {
        MASHDIST(samples)
        ch_versions = ch_versions.mix(MASHDIST.out.versions)
    } else if (params.wf == 'mashtree') {
        samples.collect{meta, fna -> fna}.map{ fna -> [[id: 'mashtree'], fna]}.set{ ch_merge_fna }
        MASHTREE(ch_merge_fna)
        ch_versions = ch_versions.mix(MASHTREE.out.versions)
    } else if (params.wf == 'mcroni') {
        MCRONI(samples)
        ch_versions = ch_versions.mix(MCRONI.out.versions)
    } else if (params.wf == 'meningotype') {
        MENINGOTYPE(samples)
        ch_versions = ch_versions.mix(MENINGOTYPE.out.versions)
    } else if (params.wf == 'merlin') {
        DATASETS()
        MERLIN(samples, DATASETS.out.mash_db)
        ch_versions = ch_versions.mix(MERLIN.out.versions)
    } else if (params.wf == 'midas') {
        MIDAS(samples)
        ch_versions = ch_versions.mix(MIDAS.out.versions)
    } else if (params.wf == 'mlst') {
        DATASETS()
        MLST(samples, DATASETS.out.mlst_db)
        ch_versions = ch_versions.mix(MLST.out.versions)
    } else if (params.wf == 'mobsuite') {
        MOBSUITE(samples)
        ch_versions = ch_versions.mix(MOBSUITE.out.versions)
    } else if (params.wf == 'mykrobe') {
        MYKROBE(samples)
        ch_versions = ch_versions.mix(MYKROBE.out.versions)
    } else if (params.wf == 'ngmaster') {
        NGMASTER(samples)
        ch_versions = ch_versions.mix(NGMASTER.out.versions)
    } else if (params.wf == 'pangenome') {
        samples.collect{meta, gff -> gff}.map{ gff -> [[id: params.use_panaroo? 'panaroo' : (params.use_roary ? 'roary' : 'pirate')], gff]}.set{ ch_merge_gff }
        PANGENOME(ch_merge_gff)
        ch_versions = ch_versions.mix(PANGENOME.out.versions)
    } else if (params.wf == 'pasty') {
        PASTY(samples)
        ch_versions = ch_versions.mix(PASTY.out.versions)
    } else if (params.wf == 'pbptyper') {
        PBPTYPER(samples)
        ch_versions = ch_versions.mix(PBPTYPER.out.versions)
    } else if (params.wf == 'plasmidfinder') {
        PLASMIDFINDER(samples)
        ch_versions = ch_versions.mix(PLASMIDFINDER.out.versions)
    } else if (params.wf == 'rgi') {
        RGI(samples)
        ch_versions = ch_versions.mix(RGI.out.versions)
    } else if (params.wf == 'seqsero2') {
        SEQSERO2(samples)
        ch_versions = ch_versions.mix(SEQSERO2.out.versions)
    } else if (params.wf == 'seroba') {
        SEROBA(samples)
        ch_versions = ch_versions.mix(SEROBA.out.versions)
    } else if (params.wf == 'shigatyper') {
        SHIGATYPER(samples)
        ch_versions = ch_versions.mix(SHIGATYPER.out.versions)
    } else if (params.wf == 'shigeifinder') {
        SHIGEIFINDER(samples)
        ch_versions = ch_versions.mix(SHIGEIFINDER.out.versions)
    } else if (params.wf == 'sistr') {
        SISTR(samples)
        ch_versions = ch_versions.mix(SISTR.out.versions)
    } else if (params.wf == 'snippy') {
        SNIPPY(samples)
        ch_versions = ch_versions.mix(SNIPPY.out.versions)
    } else if (params.wf == 'spatyper') {
        SPATYPER(samples)
        ch_versions = ch_versions.mix(SPATYPER.out.versions)
    } else if (params.wf == 'ssuissero') {
        SSUISSERO(samples)
        ch_versions = ch_versions.mix(SSUISSERO.out.versions)
    } else if (params.wf == 'staphtyper') {
        STAPHTYPER(samples)
        ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    } else if (params.wf == 'staphopiasccmec') {
        STAPHOPIASCCMEC(samples)
        ch_versions = ch_versions.mix(STAPHOPIASCCMEC.out.versions)
    } else if (params.wf == 'stecfinder') {
        STECFINDER(samples)
        ch_versions = ch_versions.mix(STECFINDER.out.versions)
    } else if (params.wf == 'tbprofiler') {
        TBPROFILER(samples)
        ch_versions = ch_versions.mix(TBPROFILER.out.versions)
    } else if (params.wf == 'tblastn') {
        TBLASTN(samples)
        ch_versions = ch_versions.mix(TBLASTN.out.versions)
    } else if (params.wf == 'tblastx') {
        TBLASTX(samples)
        ch_versions = ch_versions.mix(TBLASTX.out.versions)
    } else if (params.wf == 'teton') {
        TETON(samples)
        ch_versions = ch_versions.mix(TETON.out.versions)
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
