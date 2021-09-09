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
    if (line.runtype == 'single-end') {
        return tuple(line.sample, line.runtype, genome_size, [file(line.r1)], [params.empty_r2], file(params.empty_extra))
    } else if (line.runtype == 'paired-end') {
        return tuple(line.sample, line.runtype, genome_size, [file(line.r1)], [file(line.r2)], file(params.empty_extra))
    } else if (line.runtype == 'hybrid') {
        return tuple(line.sample, line.runtype, genome_size, [file(line.r1)], [file(line.r2)], file(line.extra))
    } else if (line.runtype == 'assembly') {
        return tuple(line.sample, line.runtype, genome_size, [params.empty_r1], [params.empty_r2], file(line.extra))
    } else if (line.runtype == 'merge-pe') {
        return tuple(line.sample, line.runtype, genome_size, handle_multiple_fqs(line.r1), handle_multiple_fqs(line.r2), file(params.empty_extra))
    } else if (line.runtype == 'hybrid-merge-pe') {
        return tuple(line.sample, line.runtype, genome_size, handle_multiple_fqs(line.r1), handle_multiple_fqs(line.r2), file(line.extra))
    } else if (line.runtype == 'merge-se') {
        return tuple(line.sample, line.runtype, genome_size, handle_multiple_fqs(line.r1), [params.empty_r2], file(params.empty_extra))
    } else {
        log.error("Invalid run_type ${line.runtype} found, please correct to continue. Expected: single-end, paired-end, hybrid, merge-pe, hybrid-merge-pe, merge-se, or assembly")
        exit 1
    }
}

def process_accessions(accession, genome_size) {
    /* Parse line and determine if single end or paired reads*/
    if (accession.length() > 0) {
        if (accession.startsWith('GCF') || accession.startsWith('GCA')) {
            return tuple(accession.split(/\./)[0], "assembly_accession", genome_size, [params.empty_r1], [params.empty_r2], file(params.empty_extra))
        } else if (accession.startsWith('DRX') || accession.startsWith('ERX') || accession.startsWith('SRX')) {
            return tuple(accession, "sra_accession", genome_size, [params.empty_r1], [params.empty_r2], file(params.empty_extra))
        } else {
            log.error("Invalid accession: ${accession} is not an accepted accession type. Accessions must be Assembly (GCF_*, GCA*) or Exeriment (DRX*, ERX*, SRX*) accessions. Please correct to continue.\n\nYou can use 'bactopia search' to convert BioProject, BioSample, or Run accessions into an Experiment accession.")
            exit 1
        }
    }
}

def create_input_channel(run_type, genome_size) {
    if (run_type == "fastqs") {
        return Channel.fromfile( file(params.fastqs) )
            .splitCsv(header: true, sep: '\t')
            .map { row -> process_fastqs(row, genome_size) }
    } else if (run_type == "is_accessions") {
        return Channel.fromfile( file(params.accessions) )
            .splitText()
            .map { line -> process_accessions(line.trim(), genome_size) }
    } else if (run_type == "is_accession") {
        return Channel.fromList([process_accessions(params.accession, log)])
    } else if (run_type == "paired-end") {
        return Channel.fromList([tuple(params.sample, run_type, genome_size, [file(params.R1)], [file(params.R2)], file(params.empty_extra))])
    } else if (run_type == "hybrid") {
        return Channel.fromList([tuple(params.sample, run_type, genome_size, [file(params.R1)], [file(params.R2)], file(params.SE))])
    } else if (run_type == "assembly") {
        return Channel.fromList([tuple(params.sample, run_type, genome_size, [params.empty_r1], [params.empty_r2], file(params.assembly))])
    } else {
        return Channel.fromList([tuple(params.sample, run_type, genome_size, [file(params.SE)], [params.empty_r2], file(params.empty_extra))])
    }
}
