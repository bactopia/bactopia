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
        def ArrayList missing_required = []
        def Integer missing_file = 0

        
        if (params.bactopia) {
            if (Utils.isLocal(params.bactopia)) {
                error += Utils.fileNotFound(params.bactopia, 'bactopia', log)
                if (error > 0) {
                    missing_required += "--bactopia"
                }
            }
        } else {
            missing_required += "--bactopia"
        }

        if (params.include && params.exclude) {
            log.error "'--include' and '--exclude' cannot be used together"
            error += 1
        } else if (params.include) {
            if (Utils.isLocal(params.include)) {
                error += Utils.fileNotFound(params.include, 'include', log)
            }
        } else if (params.exclude) {
            if (Utils.isLocal(params.exclude)) {
                error += Utils.fileNotFound(params.exclude, 'exclude', log)
            }
        }

        // Workflow specific databases
        if (params.wf == "ariba") {
            if (!params.ariba_db) {
                error += 1
                missing_required += "--ariba_db"
            }
            if (!params.ariba_dir) {
                error += 1
                missing_required += "--ariba_dir"
            }
        } else if (params.wf == "bakta") {
            if (params.bakta_db) {
                if (Utils.isLocal(params.bakta_db)) {
                    error += Utils.fileNotFound(params.bakta_db, 'bakta_db', log)
                }
            } else {
                missing_required += "--bakta_db"
            }
        } else if (params.wf == "eggnog") {
            if (params.eggnog) {
                if (Utils.isLocal(params.eggnog)) {
                    missing_file += Utils.fileNotFound("${params.eggnog}/eggnog.db", 'eggnog', log)
                    if (missing_file > 0 && params.download_eggnog == false) {
                        missing_required += "--eggnog"
                    }
                }
            } else {
                missing_required += "--eggnog"
            }
        } else if (params.wf == "gtdb") {
            if (params.gtdb) {
                if (!params.download_gtdb) {
                    if (Utils.isLocal(params.gtdb)) {
                        error += Utils.fileNotFound(params.gtdb, 'gtdb', log)
                    }
                }
            } else {
                missing_required += "--gtdb"
            }
        } else if (params.wf == "mashdist" || params.wf == "merlin") {
            if (params.mash_sketch) {
                if (Utils.isLocal(params.mash_sketch)) {
                    error += Utils.fileNotFound(params.mash_sketch, 'mash_sketch', log)
                }
            } else {
                missing_required += "--mash_sketch"
            }
        } else if (params.wf == "mykrobe") {
            if (!params.mykrobe_species) {
                error += 1
                missing_required += "--mykrobe_species"
            }
        } else if (params.wf == "pangenome") {
            if (params.traits) {
                if (Utils.isLocal(params.traits)) {
                    error += Utils.fileNotFound(params.traits, 'traits', log)
                }
            }
        } else if (params.wf == "scoary") {
            if (params.traits) {
                if (Utils.isLocal(params.traits)) {
                    error += Utils.fileNotFound(params.traits, 'traits', log)
                }
            } else {
                missing_required += "--traits"
            }
        }

        // Check for existing output directory
        if (Utils.isLocal(params.outdir)) {
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
        }

        if (missing_required.size() > 0) {
            log.error "Required parameters are missing, please check: " + missing_required.join(", ")
            log.info NfcoreSchema.paramsRequired(workflow, params, schema_filename=schema_filename)
            error += 1
        }

        if (error > 0) {
            log.error("ERROR: Validation of pipeline parameters failed!\nPlease correct to continue")
            System.exit(1)
        }
    }
}
