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
        } else if (params.wf == "bakta") {
            if (params.bakta_db) {
                if (!params.download_bakta) {
                    if (Utils.isLocal(params.bakta_db)) {
                        if (params.bakta_db.endsWith(".tar.gz")) {
                            error += Utils.fileNotFound(params.bakta_db, 'bakta_db', log)
                        } else {
                            error += Utils.fileNotFound("${params.bakta_db}/bakta.db", 'bakta_db', log)
                        }
                    }
                }
            } else {
                missing_required += "--bakta_db"
            }
        } else if (params.wf == "blastn") {
            if (params.blastn_query) {
                if (Utils.isLocal(params.blastn_query)) {
                    error += Utils.fileNotFound(params.blastn_query, 'blastn_query', log)
                }
            } else {
                missing_required += "--blastn_query"
            }
        } else if (params.wf == "blastp") {
            if (params.blastp_query) {
                if (Utils.isLocal(params.blastp_query)) {
                    error += Utils.fileNotFound(params.blastp_query, 'blastp_query', log)
                }
            } else {
                missing_required += "--blastp_query"
            }
        } else if (params.wf == "blastx") {
            if (params.blastx_query) {
                if (Utils.isLocal(params.blastx_query)) {
                    error += Utils.fileNotFound(params.blastx_query, 'blastx_query', log)
                }
            } else {
                missing_required += "--blastx_query"
            }
        } else if (params.wf == "eggnog") {
            if (params.eggnog_db) {
                if (!params.download_eggnog) {
                    if (Utils.isLocal(params.eggnog_db)) {
                        if (params.eggnog_db.endsWith(".tar.gz")) {
                            error += Utils.fileNotFound(params.eggnog_db, 'eggnog_db', log)
                        } else {
                            error += Utils.fileNotFound("${params.eggnog_db}/eggnog.db", 'eggnog_db', log)
                        }
                    }
                }
            } else {
                missing_required += "--eggnog_db"
            }
        } else if (params.wf == "gtdb") {
            if (params.gtdb) {
                if (!params.download_gtdb) {
                    if (Utils.isLocal(params.gtdb)) {
                        if (params.gtdb.endsWith(".tar.gz")) {
                            error += Utils.fileNotFound(params.gtdb, 'gtdb', log)
                        } else {
                            error += Utils.fileNotFound("${params.gtdb}/metadata/metadata.txt", 'gtdb', log)
                        }
                    }
                }
            } else {
                missing_required += "--gtdb"
            }
        } else if (params.wf == "kraken2") {
            if (params.kraken2_db) {
                if (Utils.isLocal(params.kraken2_db)) {
                    if (params.kraken2_db.endsWith(".tar.gz")) {
                        error += Utils.fileNotFound(params.kraken2_db, 'kraken2_db', log)
                    } else {
                        error += Utils.fileNotFound("${params.kraken2_db}/hash.k2d", 'kraken2_db', log)
                    }
                }
            } else {
                missing_required += "--kraken2_db"
            }
        } else if (params.wf == "mashdist") {
            if (params.mash_sketch) {
                if (Utils.isLocal(params.mash_sketch)) {
                    error += Utils.fileNotFound(params.mash_sketch, 'mash_sketch', log)
                }
            } else {
                missing_required += "--mash_sketch"
            }
        } else if (params.wf == "midas") {
            if (params.midas_db) {
                if (Utils.isLocal(params.midas_db)) {
                    if (params.midas_db.endsWith(".tar.gz")) {
                        error += Utils.fileNotFound(params.midas_db, 'midas_db', log)
                    } else {
                        error += Utils.fileNotFound("${params.midas_db}/genome_info.txt", 'midas_db', log)
                    }
                }
            } else {
                missing_required += "--midas_db"
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
        } else if (params.wf == "snippy") {
            if (params.accession && params.reference) {
                log.error "'--accession' and '--reference' cannot be used together"
                error += 1
            } else if (params.reference) {
                if (Utils.isLocal(params.reference)) {
                    error += Utils.fileNotFound(params.reference, 'reference', log)
                }
            } else if (!params.accession && !params.reference) {
                log.error "Either '--accession' and '--reference' is required"
                error += 1
            }
        } else if (params.wf == "srahumanscrubber") {
            if (params.scrubber_db) {
                if (!params.download_scrubber) {
                    if (Utils.isLocal(params.scrubber_db)) {
                        error += Utils.fileNotFound(params.scrubber_db, 'scrubber_db', log)
                    }
                }
            } else {
                missing_required += "--scrubber_db"
            }
        } else if (params.wf == "tblastn") {
            if (params.tblastn_query) {
                if (Utils.isLocal(params.tblastn_query)) {
                    error += Utils.fileNotFound(params.tblastn_query, 'tblastn_query', log)
                }
            } else {
                missing_required += "--tblastn_query"
            }
        } else if (params.wf == "tblastx") {
            if (params.tblastx_query) {
                if (Utils.isLocal(params.tblastx_query)) {
                    error += Utils.fileNotFound(params.tblastx_query, 'tblastx_query', log)
                }
            } else {
                missing_required += "--tblastx_query"
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
