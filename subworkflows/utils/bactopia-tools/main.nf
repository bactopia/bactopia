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
 * @output blastdb                BLAST database: record(meta, blastdb)
 */
nextflow.enable.types = true

include { bactopiaToolInputs } from 'plugin/nf-bactopia'
include { validateParameters } from 'plugin/nf-bactopia'

workflow BACTOPIATOOL_INIT {

    main:
    // Handle parameters
    println("Validating parameters for Bactopia Tools...")
    def validation = validateParameters(true)
    println("Validation complete.")
    if (validation.hasErrors) {
        log.info(validation.error)
        error(" ")
    } else {
        log.info(validation.logs)
    }

    // Process inputs
    def collectedInputs = bactopiaToolInputs()
    if (collectedInputs.hasErrors) {
        log.info(collectedInputs.error)
        error(" ")
    } else {
        log.info(collectedInputs.logs)
        sleep(5000)
    }

    // Build lists in a single pass (O(n) for 10k+ samples)
    def samples                    = collectedInputs.samples
    def reads_list                 = []
    def assembly_list              = []
    def assembly_reads_list        = []
    def assembly_meta_list         = []
    def assembly_proteins_gff_list = []
    def blastdb_list               = []
    def proteins_list              = []
    def gff_list                   = []
    def gbff_list                  = []

    samples.each { sample ->
    log.info("Processing sample: ${sample}")
        reads_list.add(
            record(meta: sample.meta, r1: sample.r1, r2: sample.r2, se: sample.se, lr: sample.lr)
        )
        assembly_list.add(
            record(meta: sample.meta, fna: sample.fna)
        )
        assembly_reads_list.add(
            record(meta: sample.meta, fna: sample.fna, r1: sample.r1, r2: sample.r2, se: sample.se, lr: sample.lr)
        )
        assembly_meta_list.add(
            record(meta: sample.meta, fna: sample.fna, tsv_meta: sample.tsv_meta)
        )
        assembly_proteins_gff_list.add(
            record(meta: sample.meta, fna: sample.fna_anno, faa: sample.faa, gff: sample.gff)
        )
        blastdb_list.add(
            record(meta: sample.meta, blastdb: sample.blastdb)
        )
        proteins_list.add(
            record(meta: sample.meta, faa: sample.faa)
        )
        gff_list.add(
            record(meta: sample.meta, gff: sample.gff)
        )
        gbff_list.add(
            record(meta: sample.meta, gbff: sample.gbk)
        )
    }

    emit:
    reads: Channel<Record>                 = channel.fromList(reads_list)
    assembly: Channel<Record>              = channel.fromList(assembly_list)
    assembly_reads: Channel<Record>        = channel.fromList(assembly_reads_list)
    assembly_meta: Channel<Record>         = channel.fromList(assembly_meta_list)
    assembly_proteins_gff: Channel<Record> = channel.fromList(assembly_proteins_gff_list)
    blastdb: Channel<Record>               = channel.fromList(blastdb_list)
    proteins: Channel<Record>              = channel.fromList(proteins_list)
    gff: Channel<Record>                   = channel.fromList(gff_list)
    gbff: Channel<Record>                  = channel.fromList(gbff_list)
}
