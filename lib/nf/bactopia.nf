import groovy.json.JsonSlurper

/*
========================================================================================
    Secondary input validation
========================================================================================
*/
def check_input_fastqs() {
    /* Read through --fastqs and verify each input exists. */
    USING_MERGE = false
    samples = [:]
    error = false
    has_valid_header = false
    line = 1
    file(params.fastqs).splitEachLine('\t') { cols ->
        if (line == 1) {
            if (cols[0] == 'sample' && cols[1] == 'runtype' && cols[2] == 'r1' && cols[3] == 'r2' && cols[4] == 'extra') {
                has_valid_header = true
            }
        } else {
            if (cols[0]) {
                // Sample column has a value
                if (samples.containsKey(cols[0])) {
                    samples[cols[0]] = samples[cols[0]] + 1
                } else {
                    samples[cols[0]] = 1
                }
            } else {
                log.error "LINE " + line + ':ERROR: Please verify sample name is not null'
                error = true
            }

            if (cols[2]) {
                count = 0
                cols[2].split(',').each{ fq ->
                    if (!file(fq).exists()) {
                        log.error "LINE " + line + ':ERROR: Please verify ' + fq + ' exists, and try again'
                        error = true
                    }
                    count = count + 1
                }
                if (count > 1) { 
                    USING_MERGE = true
                }
            }
            if (cols[3]) {
                cols[3].split(',').each{ fq ->
                    if (!file(fq).exists()) {
                        log.error "LINE " + line + ':ERROR: Please verify ' + fq + ' exists, and try again'
                        error = true
                    }
                }
            }
            if (cols[4]) {
                if (!file(cols[4]).exists()) {
                    log.error "LINE " + line + ':ERROR: Please verify ' + cols[4]+ ' exists, and try again'
                    error = true
                }
            }
        }
        line = line + 1
    }

    samples.each{ sample, count ->
        if (count > 1) {
            error = true
            log.error 'Sample name "'+ sample +'" is not unique, please revise sample names'
        }
    }

    if (!has_valid_header) {
        error = true
        log.error 'The header line (line 1) does not follow expected structure. (sample, runtype, r1, r2, extra)'
    }

    if (error) {
        log.error 'Verify sample names are unique and/or FASTA/FASTQ paths are correct'
        log.error 'See "--example_fastqs" for an example'
        log.error 'Exiting'
        exit 1
    }

    if (USING_MERGE) {
        log.info "\n"
        log.warn "One or more samples consists of multiple read sets that will be merged. "
        log.warn "This is an experimental feature, please use with caution."
        log.info "\n"
    }
}

def handle_multiple_fqs(read_set) {
    def fqs = []
    def String[] reads = read_set.split(",");
    reads.each { fq ->
        fqs << file(fq)
    }
    return fqs
}

def process_fastqs(line, genome_size) {
    /* Parse line and determine if single end or paired reads*/
    def meta = [:]
    meta.id = line.sample
    meta.runtype = line.runtype
    meta.genome_size = genome_size
    if (line.sample) {
        if (line.runtype == 'single-end' || line.runtype == 'ont') {
            return tuple(meta, [file(line.r1)], [params.empty_r2], file(params.empty_extra))
        } else if (line.runtype == 'paired-end') {
            return tuple(meta, [file(line.r1)], [file(line.r2)], file(params.empty_extra))
        } else if (line.runtype == 'hybrid') {
            return tuple(meta, [file(line.r1)], [file(line.r2)], file(line.extra))
        } else if (line.runtype == 'assembly') {
            return tuple(meta, [params.empty_r1], [params.empty_r2], file(line.extra))
        } else if (line.runtype == 'merge-pe') {
            return tuple(meta, handle_multiple_fqs(line.r1), handle_multiple_fqs(line.r2), file(params.empty_extra))
        } else if (line.runtype == 'hybrid-merge-pe') {
            return tuple(meta, handle_multiple_fqs(line.r1), handle_multiple_fqs(line.r2), file(line.extra))
        } else if (line.runtype == 'merge-se') {
            return tuple(meta, handle_multiple_fqs(line.r1), [params.empty_r2], file(params.empty_extra))
        } else {
            log.error("Invalid runtype ${line.runtype} found, please correct to continue. Expected: single-end, paired-end, hybrid, merge-pe, hybrid-merge-pe, merge-se, or assembly")
            exit 1
        }
    } else {
        log.error("Sample name cannot be null: ${line}")
        exit 1
    }
}

def process_accessions(accession, genome_size) {
    /* Parse line and determine if single end or paired reads*/
    def meta = [:]
    meta.genome_size = genome_size
    if (accession.length() > 0) {
        if (accession.startsWith('GCF') || accession.startsWith('GCA')) {
            meta.id = accession.split(/\./)[0]
            meta.runtype = "assembly_accession"
            return tuple(meta, [params.empty_r1], [params.empty_r2], file(params.empty_extra))
        } else if (accession.startsWith('DRX') || accession.startsWith('ERX') || accession.startsWith('SRX')) {
            meta.id = accession
            meta.runtype = "sra_accession"
            return tuple(meta, [params.empty_r1], [params.empty_r2], file(params.empty_extra))
        } else {
            log.error("Invalid accession: ${accession} is not an accepted accession type. Accessions must be Assembly (GCF_*, GCA*) or Exeriment (DRX*, ERX*, SRX*) accessions. Please correct to continue.\n\nYou can use 'bactopia search' to convert BioProject, BioSample, or Run accessions into an Experiment accession.")
            exit 1
        }
    }
}

def create_input_channel(runtype, genome_size) {
    if (runtype == "fastqs") {
        return Channel.fromPath( params.fastqs )
            .splitCsv(header: true, strip: true, sep: '\t')
            .map { row -> process_fastqs(row, genome_size) }
    } else if (runtype == "is_accessions") {
        return Channel.fromPath( params.accessions )
            .splitText()
            .map { line -> process_accessions(line.trim(), genome_size) }
    } else if (runtype == "is_accession") {
        return Channel.fromList([process_accessions(params.accession, genome_size)])
    } else {
        def meta = [:]
        meta.id = params.sample
        meta.runtype = runtype
        meta.genome_size = genome_size
        if (runtype == "paired-end") {
            return Channel.fromList([tuple(meta, [file(params.R1)], [file(params.R2)], file(params.empty_extra))])
        } else if (runtype == "hybrid") {
            return Channel.fromList([tuple(meta, [file(params.R1)], [file(params.R2)], file(params.SE))])
        } else if (runtype == "assembly") {
            return Channel.fromList([tuple(meta, [params.empty_r1], [params.empty_r2], file(params.assembly))])
        } else {
            return Channel.fromList([tuple(meta, [file(params.SE)], [params.empty_r2], file(params.empty_extra))])
        }
    }
}


/*
========================================================================================
    Import Bactopia Datasets
========================================================================================
*/
def read_json(json_file) {
    slurp = new JsonSlurper()
    return slurp.parseText(file(json_file).text)
}


def format_species(species) {
    /* Format species name to accepted format. */
    name = false
    if (species) {
        name = species.toLowerCase()
        if (species.contains(" ")) {
            name = name.replace(" ", "-")
        }
    }

    return name
}

def print_dataset_info(dataset_list, dataset_info) {
    if (dataset_list.size() > 0) {
        log.info "Found ${dataset_list.size()} ${dataset_info}"
        count = 0
        try {
            dataset_list.each {
                if (count < 5) {
                    log.info "\t${it}"
                } else {
                    log.info "\t...More than 5, truncating..."
                    throw new Exception("break")
                }
                count++
            }
        } catch (Exception e) {}
    }
}

def dataset_exists(dataset_path) {
    if (file(dataset_path).exists()) {
        return true
    } else {
        log.warn "Warning: ${dataset_path} does not exist and will not be used for analysis."
    }
}

def setup_datasets() {
    species = format_species(params.species)
    datasets = [
        'amr': [],
        'ariba': [],
        'minmer': [],
        'mlst': [],
        'references': [],
        'blast': [],
        'mapping': [],
        'proteins': file(params.empty_proteins),
        'training_set': file(params.empty_tf),
        'genome_size': params.genome_size
    ]
    genome_size = ['min': 0, 'median': 0, 'mean': 0, 'max': 0]
    if (params.datasets) {
        dataset_path = params.datasets
        available_datasets = read_json("${dataset_path}/summary.json")

        // Antimicrobial Resistance Datasets
        if (params.skip_amr) {
            log.info "Found '--skip_amr', datasets for AMRFinder+ will not be used for analysis."
        } else {
            if (available_datasets.containsKey('antimicrobial-resistance')) {
                available_datasets['antimicrobial-resistance'].each {
                    if (dataset_exists("${dataset_path}/antimicrobial-resistance/${it.name}")) {
                        datasets['amr'] << file("${dataset_path}/antimicrobial-resistance/${it.name}")
                    }
                }
                print_dataset_info(datasets['amr'], "Antimicrobial resistance datasets")
            }
        }

        // Ariba Datasets
        if (available_datasets.containsKey('ariba')) {
            available_datasets['ariba'].each {
                if (dataset_exists("${dataset_path}/ariba/${it.name}")) {
                    datasets['ariba'] << file("${dataset_path}/ariba/${it.name}")
                }
            }
            print_dataset_info(datasets['ariba'] , "ARIBA datasets")
        }

        // RefSeq/GenBank Check
        if (available_datasets.containsKey('minmer')) {
            if (available_datasets['minmer'].containsKey('sketches')) {
                available_datasets['minmer']['sketches'].each {
                    if (dataset_exists("${dataset_path}/minmer/${it}")) {
                        datasets['minmer'] << file("${dataset_path}/minmer/${it}")
                    }
                }
            }
        }
        print_dataset_info(datasets['minmer'], "minmer sketches/signatures")

        if (species) {
            if (available_datasets.containsKey('species-specific')) {
                if (available_datasets['species-specific'].containsKey(species)) {
                    species_db = available_datasets['species-specific'][species]

                    // Species-specific genome size
                    if (species_db.containsKey('genome_size')) {
                        if (['min', 'median', 'mean', 'max'].contains(datasets['genome_size'])) {
                            datasets['genome_size'] = species_db['genome_size'][params.genome_size]
                        }
                    }

                    // Prokka proteins and/or Prodigal training set
                    if (species_db.containsKey('annotation')) {
                        if (species_db['annotation'].containsKey('proteins')) {
                            prokka = "${dataset_path}/${species_db['annotation']['proteins']}"
                            if (dataset_exists(prokka)) {
                                datasets['proteins'] = file(prokka)
                                log.info "Found Prokka proteins file"
                                log.info "\t${datasets['proteins']}"
                            }
                        }

                        if (species_db['annotation'].containsKey('training_set')) {
                            prodigal_tf = "${dataset_path}/${species_db['annotation']['training_set']}"
                            if (dataset_exists(prodigal_tf)) {
                                datasets['training_set'] = file(prodigal_tf)
                                log.info "Found Prodigal training file"
                                log.info "\t${datasets['training_set'] }"
                            }
                        }
                    }

                    // Species-specific RefSeq Mash sketch for auto variant calling
                    if (!params.disable_auto_variants) {
                        if (species_db.containsKey('minmer')) {
                            if (species_db['minmer'].containsKey('mash')) {
                                refseq_minmer = "${dataset_path}/${species_db['minmer']['mash']}"
                                if (dataset_exists(refseq_minmer)) {
                                    datasets['references'] << file(refseq_minmer)
                                    log.info "Found Mash Sketch of auto variant calling"
                                    log.info "\t${refseq_minmer}"
                                }
                            }
                        }
                    }

                    // MLST 
                    if (species_db.containsKey('mlst')) {
                        species_db['mlst'].each { schema, vals ->
                            vals.each { key, val ->
                                if (key != "last_updated") {
                                    if (dataset_exists("${dataset_path}/${val}")) {
                                        datasets['mlst'] << file("${dataset_path}/${val}")
                                    }
                                }
                            }
                        }
                        print_dataset_info(datasets['mlst'], "MLST datasets")
                    }

                    if (species_db.containsKey('optional')) {
                        // References to call variants against
                        if (species_db['optional'].containsKey('reference-genomes')) {
                            file("${dataset_path}/${species_db['optional']['reference-genomes']}").list().each() {
                                if (dataset_exists("${dataset_path}/${species_db['optional']['reference-genomes']}/${it}")) {
                                    datasets['references'] << file("${dataset_path}/${species_db['optional']['reference-genomes']}/${it}")
                                }
                            }
                            print_dataset_info(datasets['references'], "reference genomes")
                        }
                        
                        // Sequences for per-base mapping
                        if (species_db['optional'].containsKey('mapping-sequences')) {
                            mapping_path = "${dataset_path}/${species_db['optional']['mapping-sequences']}"
                            mapping_total = file(mapping_path).list().size()
                            if (mapping_total > 0) {
                                datasets['mapping'] << file(mapping_path)
                            }
                            print_dataset_info(datasets['mapping'], "dir(s) of FASTAs to align reads against")
                        }

                        // FASTAs to BLAST
                        if (species_db['optional'].containsKey('blast')) {
                            species_db['optional']['blast'].each() {
                                blast_path = "${dataset_path}/${it}"
                                if (dataset_exists(blast_path)) {
                                    if (file(blast_path).list().size() > 0) {
                                        datasets['blast'] << file(blast_path)
                                    }
                                }
                            }
                            print_dataset_info(datasets['blast'], "dir(s) of FASTAs to query with BLAST")
                        }
                    }
                } else {
                    log.error "Species '${params.species}' not available, please check spelling " +
                              "or use '--available_datasets' to verify the dataset has been set " +
                              "up. Exiting"
                    exit 1
                }
            } else {
                log.error "Species '${params.species}' not available, please check spelling or " +
                          "use '--available_datasets' to verify the dataset has been set up. Exiting"
                exit 1
            }
        } else {
            log.info "--species not given, skipping the following processes (analyses):"
            log.info "\tblast (genes, proteins, primers)"
            log.info "\tcall_variants"
            log.info "\tmapping_query"
            log.info "\tsequence_type"
            if (['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
                log.error "Asked for genome size '${params.genome_size}' which requires a " +
                          "species to be given. Please give a species or specify " +
                          "a valid genome size. Exiting"
                exit 1
            }
        }
    } else {
        log.info "--datasets not given, skipping the following processes (analyses):"
        log.info "\tantimicrobial resistance"
        log.info "\tariba_analysis"
        log.info "\tblast (genes, proteins, primers)"
        log.info "\tcall_variants"
        log.info "\tmapping_query"
        log.info "\tminmer_query"
        log.info "\tsequence_type"
    }

    if (datasets['genome_size'] == null){
        log.info "Genome size not available, it will be estimated."
        datasets['genome_size'] = 0
    } else if (datasets['genome_size'].toInteger() > 0) {
        log.info "Will use ${datasets['genome_size']} bp for genome size"
    } else if (datasets['genome_size'].toInteger() == 0) {
        log.info "Found ${datasets['genome_size']} bp for genome size, it will be estimated."
    }
    
    log.info "\nIf something looks wrong, now's your chance to back out (CTRL+C 3 times). "
    log.info "Sleeping for ${params.sleep_time} seconds..."
    sleep(params.sleep_time * 1000)
    log.info "--------------------------------------------------------------------"
    return datasets
}
