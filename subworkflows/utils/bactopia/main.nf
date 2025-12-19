/**
 * Initialize Bactopia parameters and input channels.
 *
 * This subworkflow acts as the entry point for Bactopia. It uses internal plugins to
 * validate command-line parameters and organize input data (FASTQs, assemblies) into
 * a standardized channel structure for downstream analysis.
 *
 * The output uses explicit positional tuple slots with Set<Path> for each read type,
 * supporting merge operations where multiple files per slot are consolidated by GATHER:
 * - `r1`: Illumina R1 files (Set<Path>, consolidated to single file by GATHER)
 * - `r2`: Illumina R2 files (Set<Path>, consolidated to single file by GATHER)
 * - `se`: Single-end Illumina files (Set<Path>, consolidated to single file by GATHER)
 * - `lr`: Long read files (ONT/PacBio) or assembly files (Set<Path>)
 *
 * @status stable
 * @keywords initialization, validation, input-parsing, parameters, workflow
 * @tags complexity:moderate input-type:parameter output-type:single features:validation,input-parsing
 * @citation bactopia
 *
 * @input none
 * This workflow is parameter-driven and does not accept input channels.
 *
 * @output samples  Pre-merge 5-slot structure: tuple(meta, r1, r2, se, lr) as Tuple<Map, Set<Path>, Set<Path>, Set<Path>, Set<Path>>
 */
nextflow.preview.types = true

include { bactopiaInputs     } from 'plugin/nf-bactopia'
include { validateParameters } from 'plugin/nf-bactopia'

workflow BACTOPIA_INIT {

    main:
    // Handle parameters
    def validation = validateParameters(false)
    if (validation.hasErrors) {
        log.info validation.error
        error(" ")
    } else {
        log.info validation.logs
    }

    // Initialize samples channel
    def ch_samples = channel.empty() as Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>, Set<Path>>>
    def collectedInputs = bactopiaInputs(validation.data)
    if (collectedInputs.hasErrors) {
        log.info collectedInputs.error
        error(" ")
    } else {
        log.info collectedInputs.logs
    }

    collectedInputs.samples.each { sample ->
        ch_samples << tuple(
            sample.meta,
            (sample.r1 ?: []).collect { fastq -> file(fastq) }.toSet(),
            (sample.r2 ?: []).collect { fastq -> file(fastq) }.toSet(),
            (sample.se ?: []).collect { fastq -> file(fastq) }.toSet(),
            (sample.lr ?: []).collect { fastq -> file(fastq) }.toSet(),
            (sample.assembly ?: []).collect { fastq -> file(fastq) }.toSet()
        )
    }

    emit:
    // Full 6-slot structure for GATHER (pre-merge with Set<Path> for multiple files)
    samples: Channel<Tuple<Map, Set<Path?>, Set<Path?>, Set<Path?>, Set<Path?>, Set<Path?>>> = ch_samples
}
