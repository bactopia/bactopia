//
// Subworkflow with functionality specific to the Bactopia Tools
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bactopiaToolInputs } from 'plugin/nf-bactopia'
include { paramsSummaryLog   } from 'plugin/nf-bactopia'
include { validateParameters } from 'plugin/nf-bactopia'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow to initialize the Bactopia Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACTOPIATOOL_INIT {

    take:
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime

    main:

    // Handle parameters
    log.info paramsSummaryLog(workflow)
    validateParameters(null, true)

    // Collect inputs, and create appropriate tuples for 'samples' channel
    ch_samples = Channel.empty()
    bactopiaToolInputs(params.bactopia, params.workflow.ext, params.include, params.exclude).each { sample ->
        def meta = sample[0]
        def inputs = []
        def extra = []
        def extra2 = []

        // Convert string inputs to files
        if (sample[1].size() > 0) {
            sample[1].each { inputs << file(it) }
        }

        if (sample[2].size() > 0) {
            sample[2].each { extra << file(it) }
        } 

        if (sample[3].size() > 0) {
            sample[3].each { extra2 << file(it) }
            
        } 

        // Create the expected tuple
        if (extra2.size() > 0) {
            ch_samples << tuple(meta, inputs, extra, extra2)
        } else if (extra.size() > 0) {
            ch_samples << tuple(meta, inputs, extra)
        } else {
            ch_samples << tuple(meta, inputs)
        }
    }

    emit:
    samples = ch_samples
}
