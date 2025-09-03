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
include { workflowSummary   } from 'plugin/nf-bactopia'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow to initialize the Bactopia Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACTOPIATOOL_INIT {

    take:
    bactopia_path
    workflow_ext
    include_path
    exclude_path

    main:
    // Handle parameters
    log.info paramsSummaryLog(workflow)
    def validation = validateParameters(null, true)
    if (validation.hasErrors) {
        log.info validation.error
        log.info workflowSummary()
        error(" ")
    } else {
        log.info validation.logs
    }

    // Collect inputs, and create appropriate tuples for 'samples' channel
    def ch_samples = Channel.empty()
    def collectedInputs = bactopiaToolInputs(bactopia_path, workflow_ext, include_path, exclude_path)
    if (collectedInputs.hasErrors) {
        log.info collectedInputs.error
        log.info workflowSummary()
        error(" ")
    } else {
        log.info collectedInputs.logs
        sleep(5000)
    }
    collectedInputs.samples.each { sample ->
        def meta = sample[0]
        def inputs = []
        def extra = []
        def extra2 = []

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

        // Create the expected tuple
        if (extra2.size() > 0) {
            ch_samples << tuple(meta, inputs, extra, extra2)
        } else if (extra.size() > 0) {
            ch_samples << tuple(meta, inputs, extra)
        } else {
            ch_samples << tuple(meta, inputs)
        }
    }.println()

    emit:
    samples = ch_samples
}
