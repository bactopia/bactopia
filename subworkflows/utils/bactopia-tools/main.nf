//
// Subworkflow with functionality specific to the Bactopia Tools
//
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

    // Collect inputs, and create appropriate tuples for 'samples' channel
    def ch_samples = channel.empty() as Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>>
    def collectedInputs = bactopiaToolInputs()
    if (collectedInputs.hasErrors) {
        log.info collectedInputs.error
        error(" ")
    } else {
        log.info collectedInputs.logs
        sleep(5000)
    }
    collectedInputs.samples.each { sample ->
        def meta  : Map<String, String> = sample[0]
        def inputs: List<Path> = []
        def extra : List<Path> = []
        def extra2: List<Path> = []

        // Convert string inputs to files
        if (sample[1].size() > 0) {
            sample[1].each { it -> inputs << file(it) }
        }

        if (sample[2].size() > 0) {
            sample[2].each { it -> extra << file(it) }
        } 

        if (sample[3].size() > 0) {
            sample[3].each { it -> extra2 << file(it) }
        } 

        // Always create 4-element tuple with empty sets for missing data
        ch_samples << tuple(meta, inputs.toSet(), extra.toSet(), extra2.toSet())
    }

    emit:
    samples: Channel<Tuple<Map, Set<Path>>> = ch_samples.map{ meta, inputs, _extra, _extra2 ->
        tuple(meta, inputs)
    }
    samples_2: Channel<Tuple<Map, Set<Path>, Set<Path>>> = ch_samples.map{ meta, inputs, extra, _extra2 ->
        tuple(meta, inputs, extra)
    }
    samples_3: Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>> = ch_samples
}
