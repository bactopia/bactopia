/**
 * Initialize Bactopia Tool parameters and input channels.
 *
 * This subworkflow acts as the entry point for independent Bactopia Tools. It uses internal
 * plugins to validate command-line parameters and organize input data into standardized
 * channel structures suitable for various tool requirements (primary inputs only, or including
 * secondary/tertiary files).
 *
 * For FASTQ-consuming tools, provides explicit positional tuple slots:
 * - reads: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords initialization, validation, input-parsing, parameters, workflow
 * @tags complexity:moderate input-type:parameter output-type:multiple features:validation,input-parsing
 * @citation bactopia
 *
 * @input none
 * This workflow is parameter-driven and does not accept input channels.
 *
 * @output reads                  FASTQ reads: tuple(meta, r1, r2, se, lr) as Tuple<Map, Path?, Path?, Path?, Path?>
 * @output assembly               Assembly file: tuple(meta, fna) as Tuple<Map, Path>
 * @output assembly_reads         Assembly + reads: tuple(meta, fna, r1, r2, se, lr) as Tuple<Map, Path, Path?, Path?, Path?, Path?>
 * @output assembly_meta          Assembly + meta file: tuple(meta, fna, meta_file) as Tuple<Map, Path, Path>
 * @output assembly_proteins_gff  Assembly + proteins + GFF: tuple(meta, fna, faa, gff) as Tuple<Map, Path, Path?, Path?>
 * @output proteins               Protein sequences: tuple(meta, faa) as Tuple<Map, Path>
 * @output gffs                   Annotation file: tuple(meta, gff) as Tuple<Map, Path>
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
    def ch_reads                 = channel.empty() as Channel<Tuple<Map, Path?, Path?, Path?, Path?>>
    def ch_assembly              = channel.empty() as Channel<Tuple<Map, Path>>
    def ch_assembly_reads        = channel.empty() as Channel<Tuple<Map, Path, Path?, Path?, Path?, Path?>>
    def ch_assembly_meta         = channel.empty() as Channel<Tuple<Map, Path, Path>>
    def ch_assembly_proteins_gff = channel.empty() as Channel<Tuple<Map, Path, Path?, Path?>>
    def ch_proteins              = channel.empty() as Channel<Tuple<Map, Path>>
    def ch_gffs                  = channel.empty() as Channel<Tuple<Map, Path>>

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
        ch_reads                 << tuple(sample.meta, sample.r1, sample.r2, sample.se, sample.lr)
        ch_assembly              << tuple(sample.meta, sample.assembly)
        ch_assembly_reads        << tuple(sample.meta, sample.assembly, sample.r1, sample.r2, sample.se, sample.lr)
        ch_assembly_meta         << tuple(sample.meta, sample.assembly, sample.meta_file)
        ch_assembly_proteins_gff << tuple(sample.meta, sample.assembly, sample.proteins, sample.gff)
        ch_proteins              << tuple(sample.meta, sample.proteins)
        ch_gffs                  << tuple(sample.meta, sample.gff)
    }

    emit:
    reads: Channel<Tuple<Map, Path?, Path?, Path?, Path?>>                  = ch_reads
    assembly: Channel<Tuple<Map, Path>>                                     = ch_assembly
    assembly_reads: Channel<Tuple<Map, Path, Path?, Path?, Path?, Path?>>   = ch_assembly_reads
    assembly_meta: Channel<Tuple<Map, Path, Path>>                          = ch_assembly_meta
    assembly_proteins_gff: Channel<Tuple<Map, Path, Path?, Path?>>          = ch_assembly_proteins_gff
    proteins: Channel<Tuple<Map, Path>>                                     = ch_proteins
    gffs: Channel<Tuple<Map, Path>>                                         = ch_gffs
}
