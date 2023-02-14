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
        } else if  (params.R1 && params.R2 && params.SE && params.short_polish && params.sample) {
            if (Utils.isLocal(params.R1)) {
                error += Utils.fileNotGzipped(params.R1, 'R1', log)
            }
            if (Utils.isLocal(params.R2)) {
                error += Utils.fileNotGzipped(params.R2, 'R2', log)
            }
            if (Utils.isLocal(params.SE)) {
                error += Utils.fileNotGzipped(params.SE, 'SE', log)
            }
            run_type = "short_polish"
        } else if  (params.R1 && params.R2 && params.SE && params.hybrid && params.sample) {
            if (Utils.isLocal(params.R1)) {
                error += Utils.fileNotGzipped(params.R1, 'R1', log)
            }
            if (Utils.isLocal(params.R2)) {
                error += Utils.fileNotGzipped(params.R2, 'R2', log)
            }
            if (Utils.isLocal(params.SE)) {
                error += Utils.fileNotGzipped(params.SE, 'SE', log)
            }
            run_type = "hybrid"
        } else if  (params.R1 && params.R2 && params.SE) {
            log.error "Cannot use --R1, --R2, and --SE together, unless --hybrid is used."
            error += 1
        } else if  (params.R1 && params.R2 && params.sample) {
            if (Utils.isLocal(params.R1)) {
                error += Utils.fileNotGzipped(params.R1, 'R1', log)
            }
            if (Utils.isLocal(params.R2)) {
                error += Utils.fileNotGzipped(params.R2, 'R2', log)
            }
            run_type = "paired-end"
        } else if (params.SE && params.sample) {
            if (Utils.isLocal(params.SE)) {
                error += Utils.fileNotGzipped(params.SE, 'SE', log)
            }
            run_type = params.ont ? "ont" : "single-end"
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
            // The Ask Merlin feature requires a downloaded refseq mash sketch
            if (params.ask_merlin) {
                if (params.datasets) {
                    if (Utils.isLocal(params.merlin_db)) {
                        if (!Utils.fileExists(params.merlin_db)) {
                            log.error "Please verify the PATH is correct for '--merlin_db'. Unable " +
                                      "to open ${params.merlin_db}"
                            error += 1
                        }
                    }
                } else {
                    log.error "'--ask_merlin' requires '--merlin_db' to also be used"
                    error += 1
                }
            }

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

            // Check input datasets exist
            if (Utils.isLocal(params.amrfinder_db)) {
                error += Utils.fileNotFound(params.amrfinder_db, 'amrfinder_db', log)
            }
            if (Utils.isLocal(params.mlst_db)) {
                error += Utils.fileNotFound(params.mlst_db, 'mlst_db', log)
            }
        }


        // Check for existing output directory
        if (Utils.isLocal(params.outdir)) {
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
