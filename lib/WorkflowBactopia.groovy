//
// This file holds several functions specific to the workflow/bactopia.nf in Bactopia
//
// Modified from NF-Core's template: https://github.com/nf-core/tools

class WorkflowBactopia {

    //
    // Check and validate parameters
    //
    public static String initialise(workflow, params, log, schema_filename=['conf/schema/bactopia.json']) {
        def Integer error = 0
        def String run_type = ""

        if (params.fastqs) {
            error += Utils.fileNotFound(params.fastqs, 'fastqs', log)
            run_type = "fastqs"
        } else if  (params.R1 && params.R2 && params.SE && params.hybrid && params.sample) {
            error += Utils.fileNotGzipped(params.R1, 'R1', log)
            error += Utils.fileNotGzipped(params.R2, 'R2', log)
            error += Utils.fileNotGzipped(params.SE, 'SE', log)
            run_type = "hybrid"
        } else if  (params.R1 && params.R2 && params.SE) {
            log.error "Cannot use --R1, --R2, and --SE together, unless --hybrid is used."
        } else if  (params.R1 && params.R2 && params.sample) {
            error += Utils.fileNotGzipped(params.R1, 'R1', log)
            error += Utils.fileNotGzipped(params.R2, 'R2', log)
            run_type = "paired-end"
        } else if (params.SE && params.sample) {
            error += Utils.fileNotGzipped(params.SE, 'SE', log)
            run_type = params.ont ? "ont" : "single-end"
        } else if (params.assembly && params.sample) {
            error += Utils.fileNotGzipped(params.assembly, 'assembly', log)
            run_type = "assembly"
        } else if (params.accessions) {
            error += Utils.fileNotFound(params.accessions, 'accessions', log)
            run_type = "is_accessions"
        } else if (params.accession) {
            run_type = "is_accession"
        } else {
            log.error "One or more required parameters are missing, please check and try again."
            log.info NfcoreSchema.paramsRequired(workflow, params, schema_filename=schema_filename)
            error += 1
        }

        if (params.max_downloads >= 10) {
            log.warn "Please be aware the value you have set for --max_downloads (${params.max_downloads}) may cause NCBI " +
                    "to temporarily block your IP address due to too many queries at once."
        }

        if (params.genome_size) {
            if (!['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
                error += Utils.isPositiveInteger(params.genome_size, 'genome_size', log)
            }
        }

        if (params.min_time > params.max_time) {
            log.error "The value for min_time (${params.min_time}) exceeds max_time (${params.max_time}), Please correct to continue."
            error += 1
        }

        if (params.adapters) {
            error += Utils.fileNotFound(params.adapters, 'adapters', log)
        }

        if (params.phix) {
            error += Utils.fileNotFound(params.phix, 'phix', log)
        }

        if (params.datasets) {
            if (!Utils.fileExists("${params.datasets}/summary.json")) {
                log.error "Please verify the PATH is correct for '--datasets'. Unable " +
                        "to open ${params.datasets}/summary.json"
                error += 1
            }
        }

        // Check for existing output directory
        if (!params.outdir.startsWith('gs://') && !params.outdir.startsWith('s3://') && !params.outdir.startsWith('az://')) {
            // Only run this if local files
            if (!workflow.resume) {
                def run_dir = "${params.outdir}/${params.run_name}"
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
        }

        if (error > 0) {
            log.error("ERROR: Validation of pipeline parameters failed!\nPlease correct to continue")
            System.exit(1)
        }
        return run_type
    }
}
