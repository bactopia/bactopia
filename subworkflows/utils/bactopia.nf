//
// Subworkflow with functionality specific to the Bactopia Tools
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { bactopiaInputs     } from 'plugin/nf-bactopia'
include { paramsSummaryLog   } from 'plugin/nf-bactopia'
include { validateParameters } from 'plugin/nf-bactopia'
include { workflowSummary    } from 'plugin/nf-bactopia'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow to initialize the Bactopia Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACTOPIA_INIT {

    main:
    // Handle parameters
    log.info paramsSummaryLog(workflow)
    def validation = validateParameters(null, false)
    if (validation.hasErrors) {
        log.info validation.error
        log.info workflowSummary()
        error(" ")
    } else {
        log.info validation.logs
    }

    // Collect inputs, and create appropriate tuples for 'samples' channel
    def ch_samples = Channel.empty()
    def collectedInputs = bactopiaInputs(validation.data)
    if (collectedInputs.hasErrors) {
        log.info collectedInputs.error
        log.info workflowSummary()
        error(" ")
    } else {
        log.info collectedInputs.logs
    }
    log.info "${collectedInputs}"

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

        ch_samples << tuple(meta, r1, r2, extra)
    }.println()

    emit:
    samples = ch_samples
}
