import groovy.json.JsonSlurper

/*
========================================================================================
    Secondary input validation
========================================================================================
*/
def check_input_fofn() {
    /* Read through --samples and verify each input exists. */
    USING_MERGE = false
    samples = [:]
    error = false
    has_valid_header = false
    line = 1
    log.info "You provided '--check_samples', beginning to check ${params.samples}..."
    file(params.samples).splitEachLine('\t') { cols ->
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
                    if (fq.startsWith('gs:') || fq.startsWith('s3:') || fq.startsWith('az:') || fq.startsWith('https:')) {
                        log.warn "LINE " + line + ':WARN: Unable to verify existence of remote file: ' + fq
                    } else {
                        if (!file(fq).exists()) {
                            log.error "LINE " + line + ':ERROR: Please verify ' + fq + ' exists, and try again'
                            error = true
                        }
                    }

                    count = count + 1
                }
                if (count > 1) { 
                    USING_MERGE = true
                }
            }

            if (cols[3]) {
                cols[3].split(',').each{ fq ->
                    if (fq.startsWith('gs:') || fq.startsWith('s3:') || fq.startsWith('az:') || fq.startsWith('https:')) {
                        log.warn "LINE " + line + ':WARN: Unable to verify existence of remote file: ' + fq
                    } else {
                        if (!file(fq).exists()) {
                            log.error "LINE " + line + ':ERROR: Please verify ' + fq + ' exists, and try again'
                            error = true
                        }
                    }
                }
            }

            if (cols[4]) {
                    if (cols[4].startsWith('gs:') || cols[4].startsWith('s3:') || cols[4].startsWith('az:') || cols[4].startsWith('https:')) {
                        log.warn "LINE " + line + ':WARN: Unable to verify existence of remote file: ' + cols[4]
                    } else {
                        if (!file(cols[4]).exists()) {
                            log.error "LINE " + line + ':ERROR: Please verify ' + cols[4]+ ' exists, and try again'
                            error = true
                        }
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

    if (USING_MERGE) {
        log.info "\n"
        log.warn "One or more samples consists of multiple read sets that will be merged. "
        log.warn "This is an experimental feature, please use with caution."
        log.info "\n"
    }

    if (error) {
        log.error 'Verify sample names are unique and/or FASTA/FASTQ paths are correct'
        log.error 'See "bactopia prepare . --examples" for multiple example FOFNs'
        log.error 'Exiting'
        exit 1
    } else {
        log.info "Everything looked great in ${params.samples}! Feel free to start processing your genomes!"
        log.info 'Exiting'
    }
    exit 0
}

def handle_multiple_fqs(read_set) {
    def fqs = []
    def String[] reads = read_set.split(",");
    reads.each { fq ->
        fqs << file(fq)
    }
    return fqs
}

def process_fofn(line, genome_size, species) {
    /* Parse line and determine if single end or paired reads*/
    def meta = [:]
    meta.id = line.sample
    meta.runtype = line.runtype
    meta.species = line.species
    if (genome_size) {
        // User provided via --genome_size, use it
        meta.genome_size = genome_size
    } else {
        // Use size available in FOFN
        meta.genome_size = line.genome_size
    }

    if (species) {
        // User provided via --species, use it
        meta.species = species
    } else {
        // Use species available in FOFN
        meta.species = line.species
    }
    
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

def process_accessions(line, genome_size, species) {
    /* Parse line and determine if single end or paired reads*/
    def meta = [:]

    if (line.accession.startsWith('GCF') || line.accession.startsWith('GCA')) {
        meta.id = accession.split(/\./)[0]
        meta.runtype = "assembly_accession"
        meta.genome_size = genome_size
        meta.species = species
        return tuple(meta, [params.empty_r1], [params.empty_r2], file(params.empty_extra))
    } else if (line.accession.startsWith('DRX') || line.accession.startsWith('ERX') || line.accession.startsWith('SRX')) {
        meta.id = line.accession
        meta.runtype = line.runtype == 'ont' ? "sra_accession_ont" : "sra_accession"

        // If genome_size is provided, use it, otherwise use the genome_size from the FOFN
        meta.genome_size = genome_size > 0 ? genome_size : line.genome_size

        // If species is provided, use it, otherwise use the species from the FOFN
        meta.species = species ? species : line.species
    } else {
        log.error("Invalid accession: ${accession} is not an accepted accession type. Accessions must be Assembly (GCF_*, GCA*) or Exeriment (DRX*, ERX*, SRX*) accessions. Please correct to continue.\n\nYou can use 'bactopia search' to convert BioProject, BioSample, or Run accessions into an Experiment accession.")
        exit 1
    }
    return tuple(meta, [params.empty_r1], [params.empty_r2], file(params.empty_extra))
}

def process_accession(accession, genome_size, species) {
    /* Parse line and determine if single end or paired reads*/
    def meta = [:]
    meta.genome_size = genome_size
    meta.species = species
    if (accession.length() > 0) {
        if (accession.startsWith('GCF') || accession.startsWith('GCA')) {
            meta.id = accession.split(/\./)[0]
            meta.runtype = "assembly_accession"
        } else if (accession.startsWith('DRX') || accession.startsWith('ERX') || accession.startsWith('SRX')) {
            meta.id = accession
            meta.runtype = params.ont ? "sra_accession_ont" : "sra_accession"
        } else {
            log.error("Invalid accession: ${accession} is not an accepted accession type. Accessions must be Assembly (GCF_*, GCA*) or Exeriment (DRX*, ERX*, SRX*) accessions. Please correct to continue.\n\nYou can use 'bactopia search' to convert BioProject, BioSample, or Run accessions into an Experiment accession.")
            exit 1
        }
        return tuple(meta, [params.empty_r1], [params.empty_r2], file(params.empty_extra))
    }
}

def create_input_channel(runtype, genome_size, species) {
    if (runtype == "is_fofn") {
        return Channel.fromPath( params.samples )
            .splitCsv(header: true, strip: true, sep: '\t')
            .map { row -> process_fofn(row, genome_size, species) }
    } else if (runtype == "is_accessions") {
        return Channel.fromPath( params.accessions )
            .splitCsv(header: true, strip: true, sep: '\t')
            .map { row -> process_accessions(row, genome_size, species) }
    } else if (runtype == "is_accession") {
        return Channel.fromList([process_accession(params.accession, genome_size, species)])
    } else {
        def meta = [:]
        meta.id = params.sample
        meta.runtype = runtype
        meta.genome_size = genome_size
        meta.species = species
        if (runtype == "paired-end") {
            return Channel.fromList([tuple(meta, [file(params.R1)], [file(params.R2)], file(params.empty_extra))])
        } else if (runtype == "hybrid" || runtype == "short_polish") {
            return Channel.fromList([tuple(meta, [file(params.R1)], [file(params.R2)], file(params.SE))])
        } else if (runtype == "assembly") {
            return Channel.fromList([tuple(meta, [params.empty_r1], [params.empty_r2], file(params.assembly))])
        } else {
            return Channel.fromList([tuple(meta, [file(params.SE)], [params.empty_r2], file(params.empty_extra))])
        }
    }
}
