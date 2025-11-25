//
// Subworkflow with functionality specific to the Bactopia Tools
//
nextflow.preview.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bactopiaToolInputs } from 'plugin/nf-bactopia'
include { validateParameters } from 'plugin/nf-bactopia'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow to initialize the Bactopia Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACTOPIATOOL_INIT {

    main:
    // Handle parameters
    def validation = validateParameters(null, true)
    if (validation.hasErrors) {
        log.info validation.error
        error(" ")
    } else {
        log.info validation.logs
    }

    // Collect inputs, and create appropriate tuples for 'samples' channel
    def ch_samples = channel.empty() as Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>>
    def types = [false,false,false]
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
            types[0] = true
        }

        if (sample[2].size() > 0) {
            sample[2].each { it -> extra << file(it) }
            types[1] = true
        } 

        if (sample[3].size() > 0) {
            sample[3].each { it -> extra2 << file(it) }
            types[2] = true
        } 

        // Always create 4-element tuple with empty sets for missing data
        ch_samples << tuple(meta, inputs.toSet(), extra.toSet(), extra2.toSet())
    }

    emit:
    samples: Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>> = ch_samples
    data_types: Value<Integer> = types.count(true)
}
