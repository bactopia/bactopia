/**
 * Initialize Bactopia Tool parameters and input channels.
 *
 * This subworkflow acts as the entry point for independent Bactopia Tools. It uses internal
 * plugins to validate command-line parameters and organize input data into standardized
 * record-based channel structures suitable for various tool requirements.
 *
 * @status stable
 * @keywords initialization, validation, input-parsing, parameters, workflow
 * @tags complexity:moderate input-type:parameter output-type:multiple features:validation,input-parsing
 * @citation bactopia
 *
 * @input none
 * This workflow is parameter-driven and does not accept input channels.
 *
 * @output reads                  FASTQ reads: record(meta, r1, r2, se, lr)
 * @output assembly               Assembly file: record(meta, assembly)
 * @output assembly_reads         Assembly + reads: record(meta, fna, r1, r2, se, lr)
 * @output assembly_meta          Assembly + meta file: record(meta, fna, tsv_meta)
 * @output assembly_proteins_gff  Assembly + proteins + GFF: record(meta, fna, faa, gff)
 * @output proteins               Protein sequences: record(meta, proteins)
 * @output gff                    Annotation file: record(meta, gff)
 * @output gbff                   GenBank file: record(meta, gbff)
 */
nextflow.preview.types = true

include { bactopiaToolInputs } from 'plugin/nf-bactopia'
include { validateParameters } from 'plugin/nf-bactopia'

workflow BACTOPIATOOL_INIT {

    main:
    // Handle parameters
    println "Validating parameters for Bactopia Tools..."
    def validation = validateParameters(true)
    println "Validation complete."
    if (validation.hasErrors) {
        log.info validation.error
        error(" ")
    } else {
        log.info validation.logs
    }

    // Initialize channels for various output types
    def ch_reads                 = channel.empty() as Channel<Record>
    def ch_assembly              = channel.empty() as Channel<Record>
    def ch_assembly_reads        = channel.empty() as Channel<Record>
    def ch_assembly_meta         = channel.empty() as Channel<Record>
    def ch_assembly_proteins_gff = channel.empty() as Channel<Record>
    def ch_blastdb               = channel.empty() as Channel<Record>
    def ch_proteins              = channel.empty() as Channel<Record>
    def ch_gff                   = channel.empty() as Channel<Record>
    def ch_gbff                  = channel.empty() as Channel<Record>

    // Process inputs
    def collectedInputs = bactopiaToolInputs()
    if (collectedInputs.hasErrors) {
        log.info collectedInputs.error
        error(" ")
    } else {
        log.info collectedInputs.logs
        sleep(5000)
    }
    collectedInputs.samples.each { sample ->
        ch_reads                 << record(_meta: sample.meta, r1: sample.r1, r2: sample.r2, se: sample.se, lr: sample.lr)
        ch_assembly              << record(_meta: sample.meta, fna: sample.fna)
        ch_assembly_reads        << record(_meta: sample.meta, fna: sample.fna, r1: sample.r1, r2: sample.r2, se: sample.se, lr: sample.lr)
        ch_assembly_meta         << record(_meta: sample.meta, fna: sample.fna, tsv_meta: sample.tsv_meta)
        ch_assembly_proteins_gff << record(_meta: sample.meta, fna: sample.fna_anno, faa: sample.faa, gff: sample.gff)
        ch_blastdb               << record(_meta: sample.meta, blastdb: sample.blastdb)
        ch_proteins              << record(_meta: sample.meta, faa: sample.faa)
        ch_gff                   << record(_meta: sample.meta, gff: sample.gff)
        ch_gbff                  << record(_meta: sample.meta, gbff: sample.gbk)
    }

    emit:
    reads: Channel<Record>                 = ch_reads
    assembly: Channel<Record>              = ch_assembly
    assembly_reads: Channel<Record>        = ch_assembly_reads
    assembly_meta: Channel<Record>         = ch_assembly_meta
    assembly_proteins_gff: Channel<Record> = ch_assembly_proteins_gff
    blastdb: Channel<Record>               = ch_blastdb
    proteins: Channel<Record>              = ch_proteins
    gff: Channel<Record>                   = ch_gff
    gbff: Channel<Record>                  = ch_gbff
}
