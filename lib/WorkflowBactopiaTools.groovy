//
// This file holds several functions specific to the workflow/bactopiatools.nf in Bactopia
//
// Modified from NF-Core's template: https://github.com/nf-core/tools

class WorkflowBactopiaTools {

    //
    // Check and validate parameters
    //
    public static String initialise(workflow, params, log, schema_filename=['conf/schema/bactopiatools.json']) {
        def Integer error = 0

        if (params.bactopia) {
            error += Utils.fileNotFound(params.bactopia, 'bactopia', log)
        } else {
            log.error "Required parameters are missing, please check and try again."
            log.info NfcoreSchema.paramsRequired(workflow, params, schema_filename=schema_filename)
            error += 1
        }

        if (params.include && params.exclude) {
            log.error "'--include' and '--exclude' cannot be used together"
            error += 1
        } else if (params.include) {
            error += Utils.fileNotFound(params.include, 'include', log)
        } else if (params.exclude) {
            error += Utils.fileNotFound(params.exclude, 'exclude', log)
        }

        // Check file exists for certain parameters
        def Map check_file = [
            'accessions': params.accessions,
            'assemblies': params.assemblies,
            'traits': params.traits
        ]
        for ( f in check_file ) {
            if (f.value) {
                error += Utils.fileNotFound(f.value, f.key, log)
            }
        }

        // Check for existing output directory
        if (!workflow.resume) {
            def run_dir = "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}"
            def Integer files_found = 0
            new File(run_dir).eachFile { item ->
                if (item.getName() != "nf-reports") {
                    files_found += 1
                }
            }
            if (files_found > 0 && !params.force) {
                log.error("Output directory (${run_dir}) exists, ${params.wf} will not continue unless '--force' is used or a different run name (--run_name) is used.")
                error += 1
            }
        }

        if (error > 0) {
            log.error("ERROR: Validation of pipeline parameters failed!\nPlease correct to continue")
            System.exit(1)
        }
    }
}
