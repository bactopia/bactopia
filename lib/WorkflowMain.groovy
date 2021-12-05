//
// This file holds several functions specific to the main.nf workflow in Bactopia
//
// Modified from NF-Core's template: https://github.com/nf-core/tools

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            "* Bactopia\n" +
            "  https://doi.org/10.1128/mSystems.00190-20\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://bactopia.github.io/acknowledgements/"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log, schema_filename) {
        Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
        def help_string = ''
        def num_hidden = 0
        def logo_name = "bactopia"
        def command = "${workflow.manifest.name} --fastqs samples.txt --datasets datasets/ --species 'Staphylococcus aureus' -profile singularity"
        if (params.wf != "bactopia") {
            if (params.wf == "staphopia") {
                logo_name = "staphopia"
                command = "staphopia --fastqs samples.txt --datasets datasets/ --species 'Staphylococcus aureus' -profile singularity"
            } else {
                logo_name = "tools"
                command = "${workflow.manifest.name} tools ${params.wf} --bactopia /path/to/bactopia/results -profile singularity"
            }
        }
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs, logo_name, params.wf, params.workflows[params.wf].description)
        def print_example = true
        Map parsed_help = NfcoreSchema.paramsHelp(workflow, params, command, schema_filename, print_example)
        help_string += parsed_help['output']
        num_hidden += parsed_help['num_hidden']
        if (num_hidden > 0){
            help_string += colors.dim + "!! Hiding $num_hidden params, use --show_hidden_params (or --help_all) to show them !!\n" + colors.reset
        }
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Print available workflows
    //
    public static String workflows(workflow, params) {
        Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
        def wf_string = ''
        def num_hidden = 0
        def logo_name = "bactopia"
        def command = ""
        wf_string += NfcoreTemplate.logo(workflow, params.monochrome_logs, logo_name, params.wf, params.workflows[params.wf].description)
        wf_string += NfcoreSchema.listWorkflows(workflow, params)
        wf_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        wf_string += '\n' + citation(workflow) + '\n'
        wf_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return wf_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log, schema_filename) {
        def summary_log = ''
        def logo_name = "bactopia"
        if (params.wf != "bactopia") {
            if (params.wf == "staphopia") {
                logo_name = "staphopia"
            } else {
                logo_name = "tools"
            }
        }

        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs, logo_name, params.wf, params.workflows[params.wf].description)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params, schema_filename)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log, schema_filename=['nextflow_schema.json']) {
        // Print help to screen if required
        if (params.help || params.help_all) {
            log.info help(workflow, params, log, schema_filename)
            System.exit(0)
        } else if (params.list_wfs) {
            log.info workflows(workflow, params)
            System.exit(0)
        }

        if (params.validate_params) {
            // Validate workflow parameters via the JSON schema
            NfcoreSchema.validateParameters(workflow, params, log, schema_filename)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log, schema_filename)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            log.warn ""
            log.warn "Conda Disclaimer"
            log.warn ""
            log.warn "If you have access to Docker or Singularity, please consider"
            log.warn "running Bactopia using containers. The containers are less"
            log.warn "susceptible to Conda environment related issues (e.g. version"
            log.warn "conflicts) and errors caused by creation of conda environments"
            log.warn "in parallel (use '--max_cpus 1' to over come this error)."
            log.warn ""
            log.warn "To use containers, you can use the profile parameter"
            log.warn "    Docker: -profile docker"
            log.warn "    Singularity: -profile singularity"
            log.warn ""
            log.info NfcoreTemplate.dashedLine(params.monochrome_logs)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check the hostnames against configured profiles
        //NfcoreTemplate.hostName(workflow, params, log)
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttribute(params, attribute) {
        def val = ''
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
    }
}
