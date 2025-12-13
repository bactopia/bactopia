/**
 * Initialize Bactopia parameters and input channels.
 *
 * This subworkflow acts as the entry point for Bactopia. It uses internal plugins to
 * validate command-line parameters and organize input data (FASTQs, assemblies) into
 * a standardized channel structure for downstream analysis.
 *
 * @status stable
 * @keywords initialization, validation, input-parsing, parameters, workflow
 * @tags complexity:moderate input-type:parameter output-type:single features:validation,input-parsing
 * @citation bactopia
 *
 * @input none
 * This workflow is parameter-driven and does not accept input channels.
 *
 * @output samples        Channel containing standardized sample inputs (meta, r1, r2, extra)
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

    // Collect inputs, and create appropriate tuples for 'samples' channel
    def ch_samples = channel.empty() as Channel<Tuple<Map, Set<Path>, Set<Path>, Path>>
    def collectedInputs = bactopiaInputs(validation.data)
    if (collectedInputs.hasErrors) {
        log.info collectedInputs.error
        error(" ")
    } else {
        log.info collectedInputs.logs
    }

    collectedInputs.samples.each { sample ->
        def meta = sample[0]
        def r1 = []
        def r2 = []
        def extra = file(sample[3])

        // Convert string inputs to files
        if (sample[1].size() > 0) {
            sample[1].each { it -> r1 << file(it) }
        }

        if (sample[2].size() > 0) {
            sample[2].each { it -> r2 << file(it) }
        } 

        ch_samples << tuple(meta, r1.toSet(), r2.toSet(), extra)
    }.println()

    emit:
    samples: Channel<Tuple<Map, Set<Path>, Set<Path>, Path>> = ch_samples
}
