/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * This subworkflow orchestrates the execution of abricate components.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation abricate
 *
 * @subworkflows bactopiainputs, validateparameters
 * @output samples Samples
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
