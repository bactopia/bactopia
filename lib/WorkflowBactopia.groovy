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

        if (params.samples) {
            if (Utils.isLocal(params.samples)) {
                error += Utils.fileNotFound(params.samples, 'samples', log)
            }
            run_type = "is_fofn"
        } else if  (params.r1 && params.r2 && params.ont && params.short_polish && params.sample) {
            if (Utils.isLocal(params.r1)) {
                error += Utils.fileNotGzipped(params.r1, 'r1', log)
            }
            if (Utils.isLocal(params.r1)) {
                error += Utils.fileNotGzipped(params.r2, 'r2', log)
            }
            if (Utils.isLocal(params.ont)) {
                error += Utils.fileNotGzipped(params.ont, 'ont', log)
            }
            run_type = "short_polish"
        } else if  (params.r1 && params.r2 && params.ont && params.hybrid && params.sample) {
            if (Utils.isLocal(params.r1)) {
                error += Utils.fileNotGzipped(params.r1, 'r1', log)
            }
            if (Utils.isLocal(params.r2)) {
                error += Utils.fileNotGzipped(params.r2, 'r2', log)
            }
            if (Utils.isLocal(params.ont)) {
                error += Utils.fileNotGzipped(params.ont, 'ont', log)
            }
            run_type = "hybrid"
        } else if  (params.r1 && params.r2 && params.se) {
            log.error "Cannot use --r1, --r2, and --se together"
            error += 1
        } else if  (params.r1 && params.r2 && params.ont) {
            log.error "Cannot use --r1, --r2, and --ont together, unless using --short_polish or --hybrid"
            error += 1
        } else if  (params.ont && params.se) {
            log.error "Cannot use --ont and --se together"
            error += 1
        }  else if  (params.r1 && params.r2 && params.sample) {
            if (Utils.isLocal(params.r1)) {
                error += Utils.fileNotGzipped(params.r1, 'r1', log)
            }
            if (Utils.isLocal(params.r2)) {
                error += Utils.fileNotGzipped(params.r2, 'r2', log)
            }
            run_type = "paired-end"
        } else if (params.ont && params.sample) {
            if (Utils.isLocal(params.ont)) {
                error += Utils.fileNotGzipped(params.ont, 'ont', log)
            }
            run_type = "ont" 
        } else if (params.se && params.sample) {
            if (Utils.isLocal(params.se)) {
                error += Utils.fileNotGzipped(params.se, 'se', log)
            }
            run_type = "single-end"
        } else if (params.assembly && params.sample) {
            if (Utils.isLocal(params.assembly)) {
                error += Utils.fileNotGzipped(params.assembly, 'assembly', log)
            }
            run_type = "assembly"
        } else if (params.accessions) {
            if (Utils.isLocal(params.accessions)) {
                error += Utils.fileNotFound(params.accessions, 'accessions', log)
            }
            run_type = "is_accessions"
        } else if (params.accession) {
            run_type = "is_accession"
        } else {
            log.error "One or more required parameters are missing, please check and try again."
            log.info NfcoreSchema.paramsRequired(workflow, params, schema_filename=schema_filename)
            error += 1
        }

        if (params.check_samples && !params.samples) {
            log.error "To use --check_samples, you must also provide a FOFN to check using --samples."
            error += 1
        }

        if (params.max_downloads >= 10) {
            log.warn "Please be aware the value you have set for --max_downloads (${params.max_downloads}) may cause NCBI " +
                     "to temporarily block your IP address due to too many queries at once."
        }

        if (params.genome_size) {
            error += Utils.isPositiveInteger(params.genome_size, 'genome_size', log)
        }

        if (params.min_time > params.max_time) {
            log.error "The value for min_time (${params.min_time}) exceeds max_time (${params.max_time}), Please correct to continue."
            error += 1
        }

        if (params.containsKey('adapters')){
            if (params.adapters) {
                if (Utils.isLocal(params.adapters)) {
                    error += Utils.fileNotFound(params.adapters, 'adapters', log)
                }
            }

            if (params.phix) {
                if (Utils.isLocal(params.phix)) {
                    error += Utils.fileNotFound(params.phix, 'phix', log)
                }
            }
        }

        // following should only be checked for specific workflows
        if (['bactopia', 'staphopia'].contains(params.wf)) {
            // Using Bakta, requires path to database
            if (params.use_bakta) {
                if (params.bakta_db) {
                    if (Utils.isLocal(params.bakta_db)) {
                        if (params.bakta_db.endsWith(".tar.gz")) {
                            error += Utils.fileNotFound(params.bakta_db, 'bakta_db', log)
                        } else {
                            error += Utils.fileNotFound("${params.bakta_db}/bakta.db", 'bakta_db', log)
                        }
                    }
                } else {
                    log.error "'--use_bakta' requires '--bakta_db' to also be used"
                    error += 1
                }
            }
        }

        // Check for existing output directory
        /*
        if (Utils.isLocal(params.outdir)) {
            // Only run this if local files
            if (!workflow.resume) {
                def Integer files_found = 0
                new File("${params.outdir}/bactopia-comparative/${params.wf}/${params.run_name}").eachDirRecurse { item ->
                    if (item.toString().contains("nf-reports")) {
                        return
                    } else {
                        files_found += 1
                    }
                }

                if (files_found > 0 && !params.force) {
                    log.error("Output for ${params.run_name} (--run_name) already exists in ${params.outdir} (--outdir), ${params.wf} will not continue unless '--force' is used, a different run name (--run_name), or a different output directory (--outdir) is used.")
                    error += 1
                }
            }
        }
        */

        if (error > 0) {
            log.error("ERROR: Validation of pipeline parameters failed!\nPlease correct to continue")
            System.exit(1)
        }

        return run_type
    }
}
