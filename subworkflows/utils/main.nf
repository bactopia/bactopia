//
// Subworkflow with functionality specific to the nf-core/sarek pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLESHEET_TO_CHANNEL    } from '../samplesheet_to_channel'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion        } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(nextflow_cli_args)

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // Check input path parameters to see if they exist
    def checkPathParamList = [
        params.ascat_alleles,
        params.ascat_loci,
        params.ascat_loci_gc,
        params.ascat_loci_rt,
        params.bwa,
        params.bwamem2,
        params.bcftools_annotations,
        params.bcftools_annotations_tbi,
        params.bcftools_header_lines,
        params.cf_chrom_len,
        params.chr_dir,
        params.cnvkit_reference,
        params.dbnsfp,
        params.dbnsfp_tbi,
        params.dbsnp,
        params.dbsnp_tbi,
        params.dict,
        params.dragmap,
        params.fasta,
        params.fasta_fai,
        params.germline_resource,
        params.germline_resource_tbi,
        params.input,
        params.intervals,
        params.known_indels,
        params.known_indels_tbi,
        params.known_snps,
        params.known_snps_tbi,
        params.mappability,
        params.multiqc_config,
        params.ngscheckmate_bed,
        params.pon,
        params.pon_tbi,
        params.sentieon_dnascope_model,
        params.spliceai_indel,
        params.spliceai_indel_tbi,
        params.spliceai_snv,
        params.spliceai_snv_tbi
    ]

// only check if we are using the tools
if (params.tools && (params.tools.split(',').contains('snpeff') || params.tools.split(',').contains('merge'))) checkPathParamList.add(params.snpeff_cache)
if (params.tools && (params.tools.split(',').contains('vep')    || params.tools.split(',').contains('merge'))) checkPathParamList.add(params.vep_cache)

    // def retrieveInput(need_input, step, outdir) {

    params.input_restart = retrieveInput((!params.build_only_index && !params.input), params.step, params.outdir)

    ch_from_samplesheet = params.build_only_index ? Channel.empty() : params.input ?
        Channel.fromList(samplesheetToList(params.input, "$projectDir/assets/schema_input.json")) :
        Channel.fromList(samplesheetToList(params.input_restart, "$projectDir/assets/schema_input.json"))

    SAMPLESHEET_TO_CHANNEL(
        ch_from_samplesheet,
        params.aligner,
        params.ascat_alleles,
        params.ascat_loci,
        params.ascat_loci_gc,
        params.ascat_loci_rt,
        params.bcftools_annotations,
        params.bcftools_annotations_tbi,
        params.bcftools_header_lines,
        params.build_only_index,
        params.dbsnp,
        params.fasta,
        params.germline_resource,
        params.intervals,
        params.joint_germline,
        params.joint_mutect2,
        params.known_indels,
        params.known_snps,
        params.no_intervals,
        params.pon,
        params.sentieon_dnascope_emit_mode,
        params.sentieon_haplotyper_emit_mode,
        params.seq_center,
        params.seq_platform,
        params.skip_tools,
        params.snpeff_cache,
        params.snpeff_db,
        params.step,
        params.tools,
        params.umi_read_structure,
        params.wes)

    emit:
    samplesheet = SAMPLESHEET_TO_CHANNEL.out.input_sample
    versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    def multiqc_report_list = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_report_list.getVal()
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Bactopia Logos
//
def bactopiaLogo(workflow, monochrome_logs=true, logo_name="bactopia", worflow_name="bactopia", worflow_description="bactopia") {
    Map colors = logColours(monochrome_logs)
    if (logo_name == "cleanyerreads") {
        String.format(
            """\n
            -${colors.dim}------------------------------------------------------------------------${colors.reset}-
            ${colors.blue}    ____ _                   __   __          ____                _              ${colors.reset}
            ${colors.blue}   / ___| | ___  __ _ _ __   \\ \\ / /__ _ __  |  _ \\ ___  __ _  __| |___       ${colors.reset}
            ${colors.blue}  | |   | |/ _ \\/ _` | '_ \\   \\ V / _ \\ '__| | |_) / _ \\/ _` |/ _` / __|    ${colors.reset}
            ${colors.blue}  | |___| |  __/ (_| | | | |   | |  __/ |    |  _ <  __/ (_| | (_| \\__ \\       ${colors.reset}
            ${colors.blue}   \\____|_|\\___|\\__,_|_| |_|   |_|\\___|_|    |_| \\_\\___|\\__,_|\\__,_|___/ ${colors.reset}
            ${colors.blue} ${colors.reset}
            ${colors.blue}                      ⠀⠀⡶⠛⠲⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⡶⠚⢲⡀⠀           ${colors.reset}
            ${colors.blue}                      ⣰⠛⠃⠀⢠⣏⠀⠀⠀⠀⣀⣠⣤⣤⣤⣤⣤⣤⣤⣀⡀⠀⠀⠀⣸⡇⠀⠈⠙⣧         ${colors.reset}
            ${colors.blue}                      ⠸⣦⣤⣄⠀⠙⢷⣤⣶⠟⠛⢉⣁⣠⣤⣤⣤⣀⣉⠙⠻⢷⣤⡾⠋⢀⣠⣤⣴⠟        ${colors.reset}
            ${colors.blue}                      ⠀⠀⠀⠈⠳⣤⡾⠋⣀⣴⣿⣿⠿⠿⠟⠛⠿⠿⣿⣿⣶⣄⠙⢿⣦⠟⠁⠀⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⠀⢀⣾⠟⢀⣼⣿⠟⠋⠀⠀⠀⠀⠀⠀⠀⠀⠉⠻⣿⣷⡄⠹⣷⡀⠀⠀⠀          ${colors.reset}
            ${colors.blue}                      ⠀⠀⠀⣾⡏⢠⣿⣿⡯⠤⠤⠤⠒⠒⠒⠒⠒⠒⠒⠤⠤⠽⣿⣿⡆⠹⣷⡀⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⢸⣟⣠⡿⠿⠟⠒⣒⣒⣈⣉⣉⣉⣉⣉⣉⣉⣁⣒⣒⡛⠻⠿⢤⣹⣇⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⣾⡭⢤⣤⣠⡞⠉⠉⢀⣀⣀⠀⠀⠀⠀⢀⣀⣀⠀⠈⢹⣦⣤⡤⠴⣿⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⣿⡇⢸⣿⣿⣇⠀⣼⣿⣿⣿⣷⠀⠀⣼⣿⣿⣿⣷⠀⢸⣿⣿⡇⠀⣿⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⢻⡇⠸⣿⣿⣿⡄⢿⣿⣿⣿⡿⠀⠀⢿⣿⣿⣿⡿⢀⣿⣿⣿⡇⢸⣿⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⠸⣿⡀⢿⣿⣿⣿⣆⠉⠛⠋⠁⢴⣶⠀⠉⠛⠉⣠⣿⣿⣿⡿⠀⣾⠇⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⠀⢻⣷⡈⢻⣿⣿⣿⣿⣶⣤⣀⣈⣁⣀⡤⣴⣿⣿⣿⣿⡿⠁⣼⠟⠀⠀⠀         ${colors.reset}
            ${colors.blue}                      ⠀⠀⠀⢀⣽⣷⣄⠙⢿⣿⣿⡟⢲⠧⡦⠼⠤⢷⢺⣿⣿⡿⠋⣠⣾⢿⣄⠀⠀⠀         ${colors.reset}
            ${colors.blue}                      ⢰⠟⠛⠟⠁⣨⡿⢷⣤⣈⠙⢿⡙⠒⠓⠒⠓⠚⣹⠛⢉⣠⣾⠿⣧⡀⠙⠋⠙⣆         ${colors.reset}
            ${colors.blue}                      ⠹⣄⡀⠀⠐⡏⠀⠀⠉⠛⠿⣶⣿⣦⣤⣤⣤⣶⣷⡾⠟⠋⠀⠀⢸⡇⠀⢠⣤⠟         ${colors.reset}
            ${colors.blue}                      ⠀⠀⠳⢤⠼⠃⠀⠀⠀⠀⠀⠀⠈⠉⠉⠉⠉⠁⠀⠀⠀⠀⠀⠀⠘⠷⢤⠾⠁⠀         ${colors.reset}
            ${colors.blue} ${colors.reset}
            ${colors.cyan} clean-yer-reads v${workflow.manifest.version}${colors.reset}
            ${colors.cyan} ${worflow_description} ${colors.reset}
            -${colors.dim}------------------------------------------------------------------------${colors.reset}-
            """.stripIndent()
        )
    } else if (logo_name == "enteropia") {
        String.format(
            """\n
            -${colors.dim}---------------------------------------------------------------------------------------${colors.reset}-
            ${colors.blue}              _                       _             ${colors.reset}
            ${colors.blue}    ___ _ __ | |_ ___ _ __ ___  _ __ (_) __ _       ${colors.reset}
            ${colors.blue}   / _ \\ '_ \\| __/ _ \\ '__/ _ \\| '_ \\| |/ _` | ${colors.reset}
            ${colors.blue}  |  __/ | | | ||  __/ | | (_) | |_) | | (_| |      ${colors.reset}
            ${colors.blue}   \\___|_| |_|\\__\\___|_|  \\___/| .__/|_|\\__,_| ${colors.reset}
            ${colors.blue}                               |_|                  ${colors.reset}
            ${colors.cyan}  enteropia v${workflow.manifest.version}${colors.reset}
            ${colors.cyan}  ${worflow_description} ${colors.reset}
            -${colors.dim}---------------------------------------------------------------------------------------${colors.reset}-
            """.stripIndent()
        )
    } else if (logo_name == "staphopia") {
        String.format(
            """\n
            -${colors.dim}------------------------------------------------${colors.reset}-
            ${colors.blue}       _              _                 _            ${colors.reset}
            ${colors.blue}   ___| |_ __ _ _ __ | |__   ___  _ __ (_) __ _      ${colors.reset}
            ${colors.blue}  / __| __/ _` | '_ \\| '_ \\ / _ \\| '_ \\| |/ _` | ${colors.reset}
            ${colors.blue}  \\__ \\ || (_| | |_) | | | | (_) | |_) | | (_| |   ${colors.reset}
            ${colors.blue}  |___/\\__\\__,_| .__/|_| |_|\\___/| .__/|_|\\__,_| ${colors.reset}
            ${colors.blue}               |_|               |_|                 ${colors.reset}
            ${colors.cyan}  staphopia v${workflow.manifest.version}${colors.reset}
            ${colors.cyan}  ${worflow_description} ${colors.reset}
            -${colors.dim}------------------------------------------------${colors.reset}-
            """.stripIndent()
        )
    } else if (logo_name == "teton") {
        String.format(
            """\n
            -${colors.dim}------------------------------------------------------------------${colors.reset}-
            ${colors.blue}   _       _                          _     *                     ${colors.reset}
            ${colors.blue}  | |_ ___| |_ ___  _ __       *     / \\_       *   /\\'__       ${colors.reset}
            ${colors.blue}  | __/ _ \\ __/ _ \\| '_ \\       /\\ _/    \\        _/  /  \\  * ${colors.reset}
            ${colors.blue}  | ||  __/ || (_) | | | |     /\\/\\  /\\/  \\_   _^/  ^/    `--.  ${colors.reset}
            ${colors.blue}   \\__\\___|\\__\\___/|_| |_|    /    \\/  \\    \\ /.' ^_   \\_   .'\\ ${colors.reset}
            ${colors.blue}                                      Art by Joan Stark         ${colors.reset}
            ${colors.blue}  ${colors.reset}
            ${colors.cyan}  teton v${workflow.manifest.version}${colors.reset}
            ${colors.cyan}  ${worflow_description}${colors.reset}
            -${colors.dim}------------------------------------------------------------------${colors.reset}-
            """.stripIndent()
        )
    } else if (logo_name == "tools") {
        String.format(
            """\n
            -${colors.dim}------------------------------------------------------------------${colors.reset}-
            ${colors.blue}   _                _              _         _              _              ${colors.reset}
            ${colors.blue}  | |__   __ _  ___| |_ ___  _ __ (_) __ _  | |_ ___   ___ | |___          ${colors.reset}
            ${colors.blue}  | '_ \\ / _` |/ __| __/ _ \\| '_ \\| |/ _` | | __/ _ \\ / _ \\| / __|    ${colors.reset}
            ${colors.blue}  | |_) | (_| | (__| || (_) | |_) | | (_| | | || (_) | (_) | \\__ \\        ${colors.reset}
            ${colors.blue}  |_.__/ \\__,_|\\___|\\__\\___/| .__/|_|\\__,_|  \\__\\___/ \\___/|_|___/ ${colors.reset}
            ${colors.blue}                            |_|                                            ${colors.reset}
            ${colors.cyan}  ${workflow.manifest.name} tools ${worflow_name} v${workflow.manifest.version}${colors.reset}
            ${colors.cyan}  ${worflow_description} ${colors.reset}
            -${colors.dim}------------------------------------------------------------------${colors.reset}-
            """.stripIndent()
        )
    } else {
        String.format(
            """\n
            -${colors.dim}-------------------------------------------${colors.reset}-
            ${colors.blue}   _                _              _             ${colors.reset}
            ${colors.blue}  | |__   __ _  ___| |_ ___  _ __ (_) __ _       ${colors.reset}
            ${colors.blue}  | '_ \\ / _` |/ __| __/ _ \\| '_ \\| |/ _` |   ${colors.reset}
            ${colors.blue}  | |_) | (_| | (__| || (_) | |_) | | (_| |      ${colors.reset}
            ${colors.blue}  |_.__/ \\__,_|\\___|\\__\\___/| .__/|_|\\__,_| ${colors.reset}
            ${colors.blue}                            |_|                  ${colors.reset}
            ${colors.cyan}  ${workflow.manifest.name} v${workflow.manifest.version}${colors.reset}
            ${colors.cyan}  ${worflow_description} ${colors.reset}
            -${colors.dim}-------------------------------------------${colors.reset}-
            """.stripIndent()
        )
    }
}



//
// retrieveInput
//
def retrieveInput(need_input, step, outdir) {
    def input = null
    if (!params.input && !params.build_only_index) {
        switch (step) {
            case 'mapping':                 error("Can't start $step step without samplesheet")
                                            break
            case 'markduplicates':          log.warn("Using file ${outdir}/csv/mapped.csv");
                                            input = outdir + "/csv/mapped.csv"
                                            break
            case 'prepare_recalibration':   log.warn("Using file ${outdir}/csv/markduplicates_no_table.csv");
                                            input = outdir + "/csv/markduplicates_no_table.csv"
                                            break
            case 'recalibrate':             log.warn("Using file ${outdir}/csv/markduplicates.csv");
                                            input = outdir + "/csv/markduplicates.csv"
                                            break
            case 'variant_calling':         log.warn("Using file ${outdir}/csv/recalibrated.csv");
                                            input = outdir + "/csv/recalibrated.csv"
                                            break
            // case 'controlfreec':         csv_file = file("${outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
            case 'annotate':                log.warn("Using file ${outdir}/csv/variantcalled.csv");
                                            input = outdir + "/csv/variantcalled.csv"
                                            break
            default:                        log.warn("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
                                            error("Unknown step $step")
        }
    }
    return input
}
