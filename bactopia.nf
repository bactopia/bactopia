#! /usr/bin/env nextflow
import groovy.json.JsonSlurper
PROGRAM_NAME = 'bactopia'
VERSION = '0.0.1'
if (params.help || params.full_usage) print_usage();
if (workflow.commandLine.endsWith(workflow.scriptName)) print_usage();
if (params.example_fastqs) print_example_fastqs();
if (params.version) print_version();
check_input_params()
check_input_fastqs(params.fastqs)
if (params.check_fastqs) print_check_fastqs(params.fastqs);
if (params.available_databases) print_available_databases(params.database)

// Set the maximum number of cpus to use
config.poolSize = params.max_cpus.toInteger()
cpus = params.cpus
if (cpus > params.max_cpus) {
    cpus = params.max_cpus
    log.info "--cpus ${params.cpus} exceeded --max_cpus ${params.max_cpus}, changed ${params.cpus} to ${cpus}"
}

// Setup output directories
outdir = params.outdir ? params.outdir : './'

// Setup some defaults
ARIBA_DATABASES = []
MINMER_DATABASES = []
MLST_DATABASES = []
REFERENCES = []
INSERTIONS = []
PRIMERS = []
PROKKA_PROTEINS = null
organism_genome_size = ['min': 0, 'median': 0, 'mean': 0, 'max': 0]
if (params.database) {
    database_path = get_absolute_path(params.database)
    available_databases = read_database_summary(database_path)

    available_databases['ariba'].each {
        ARIBA_DATABASES << "${database_path}/ariba/${it.name}"
    }
    print_database_info(ARIBA_DATABASES, "ARIBA databases")

    available_databases['minmer']['sketches'].each {
        MINMER_DATABASES << "${database_path}/minmer/${it}"
    }
    print_database_info(MINMER_DATABASES, "minmer sketches databases")

    if (params.organism) {
        if (available_databases.containsKey(params.organism)) {
            organism_genome_size = available_databases[params.organism]['genome_size']
            prokka = "${database_path}/${available_databases[params.organism]['prokka']['proteins']}"
            if (file(prokka).exists()) {
                PROKKA_PROTEINS = prokka
                log.info "Found Prokka proteins files (${PROKKA_PROTEINS})"

            }
            available_databases[params.organism]['mlst'].each { key, val ->
                if (key != "last_updated") {
                    if (file("${database_path}/${val}").exists()) {
                        MLST_DATABASES << "${database_path}/${val}"
                    }
                }
            }
            print_database_info(MLST_DATABASES, "MLST databases")

            file("${database_path}/${available_databases[params.organism]['reference-genomes']}").list().each() {
                REFERENCES << "${database_path}/${params.organism}/reference-genomes/${it}"
            }
            print_database_info(REFERENCES, "reference genomes")

            file("${database_path}/${available_databases[params.organism]['insertion-sequences']}").list().each() {
                INSERTIONS << "${database_path}/${params.organism}/insertion-sequences/${it}"
            }
            print_database_info(INSERTIONS, "insertion sequence FASTAs")

            PRIMERS = file("${database_path}/${available_databases[params.organism]['primer-sequences']}").list()
            print_database_info(PRIMERS, "primer sequence FASTAs")
        } else {
            log.info "Organism '${params.organism}' not available, please check spelling or use '--available_databases' " +
                     "to verify the database has been set up. Exiting"
            exit 1
        }
    } else {
        log.info "--organism not given, skipping the following processes (analyses):"
        log.info "\tsequence_type"
        log.info "\tcall_variants"
        log.info "\tinsertion_sequence_query"
        log.info "\tprimer_query"
        if (['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
            log.info "Asked for genome size '${params.genome_size}' which requires an " +
                     "organism to be given. Please give an organism or specify " +
                     "a valid genome size. Exiting"
            exit 1
        }
    }
} else {
    log.info "--database not given, skipping the following processes (analyses):"
    log.info "\tsequence_type"
    log.info "\tariba_databases"
    log.info "\tminmer_query"
    log.info "\tcall_variants"
    log.info "\tinsertion_sequence_query"
    log.info "\tprimer_query"
}


process estimate_genome_size {
    /* Estimate the input genome size if not given. */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from create_fastq_channel(params.fastqs)

    output:
    file "genome-size.txt" into GS_QC_READS, GS_ASSEMBLY

    shell:
    template(task.ext.template)
}


process qc_reads {
    /* Cleanup the reads using Illumina-Cleanup */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from create_fastq_channel(params.fastqs)
    file(genome_size_file) from GS_QC_READS

    output:
    file "quality-control/*"
    set val(sample), val(single_end),
        file("quality-control/${sample}*.fastq.gz") into ASSEMBLY, SEQUENCE_TYPE, COUNT_31MERS,
                                                         ARIBA_ANALYSIS, MINMER_SKETCH, MINMER_QUERY,
                                                         INSERTION_SEQUENCES, CALL_VARIANTS

    shell:
    fq2 = single_end ? "" : fq[1]
    template(task.ext.template)
}


process assemble_genome {
    /* Assemble the genome using Shovill, SKESA is used by default */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/assembly", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from ASSEMBLY
    file(genome_size_file) from GS_ASSEMBLY

    output:
    file "shovill*"
    file "${sample}.fna.gz" into SEQUENCE_TYPE_ASSEMBLY
    file "${sample}.fna.json"
    set val(sample), file("${sample}.fna.gz") into ANNOTATION

    shell:
    template(task.ext.template)
}


process annotate_genome {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), file(fasta) from ANNOTATION

    output:
    file 'annotation/*'
    file 'annotation/*.gbk.gz' into INSERTION_GENBANK

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    proteins = PROKKA_PROTEINS ? "--proteins ${PROKKA_PROTEINS}" : ""
    template(task.ext.template)
}


process count_31mers {
    /* Count 31mers in the reads using McCortex */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/kmers", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from COUNT_31MERS

    output:
    file "${sample}.ctx"

    shell:
    template(task.ext.template)

}


process sequence_type {
    /* Determine MLST types using ARIBA and BLAST */
    cpus 1
    errorStrategy 'retry'
    maxRetries 5
    tag "${sample} - ${method}"
    publishDir "${outdir}/${sample}/mlst", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from SEQUENCE_TYPE
    file(assembly) from SEQUENCE_TYPE_ASSEMBLY
    each database from MLST_DATABASES

    when:
    database.contains('blast') || (database.contains('ariba') && single_end == false)

    output:
    file "${method}/*"

    shell:
    method = database.contains('blast') ? 'blast' : 'ariba'
    template(task.ext.template)
}


process ariba_analysis {
    /* Run reads against all available (if any) ARIBA databases */
    cpus 1
    errorStrategy 'retry'
    maxRetries 5
    tag "${sample} - ${database_name}"
    publishDir "${outdir}/${sample}/ariba", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from ARIBA_ANALYSIS
    each database from ARIBA_DATABASES

    output:
    file "${database_name}/*"

    when:
    single_end == false

    shell:
    database_name = file(database).getName()
    template(task.ext.template)
}



process minmer_sketch {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31) and
    Sourmash (k=21,31,51)
    */
    cpus 1
    tag "${sample}"
    publishDir "${outdir}/${sample}/minmers", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from MINMER_SKETCH

    output:
    file("${sample}*.{msh,sig}")
    file("${sample}.sig") into QUERY_SOURMASH

    shell:
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template(task.ext.template)
}



process minmer_query {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    cpus 1
    tag "${sample} - ${minmer_database}"
    publishDir "${outdir}/${sample}/minmers", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from MINMER_QUERY
    file(sourmash) from QUERY_SOURMASH
    each database from MINMER_DATABASES

    output:
    file("${sample}*.txt")

    shell:
    minmer_database = file(database).getName()
    mash_w = params.screen_w ? "-w" : ""
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template(task.ext.template)
}


process call_variants {
    /*
    Identify variants (SNPs/InDels) against a set of reference genomes
    using Snippy.
    */
    cpus cpus
    tag "${sample} - ${reference_name}"
    publishDir "${outdir}/${sample}/variants", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from CALL_VARIANTS
    each reference from REFERENCES

    output:
    file("${reference_name}/*")

    shell:
    reference_name = file(reference).getSimpleName()
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    template(task.ext.template)
}


process insertion_sequences {
    /*
    Query a set of insertion sequences (FASTA) against annotated GenBank file
    using ISMapper.
    */
    cpus cpus
    tag "${sample} - ${insertion_name}"
    publishDir "${outdir}/${sample}/", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from INSERTION_SEQUENCES
    file(genbank) from INSERTION_GENBANK
    each insertion_fasta from INSERTIONS

    output:
    file("insertion-sequences/*")

    when:
    single_end == false

    shell:
    insertion_name = file(insertion_fasta).getSimpleName()
    gunzip_genbank = genbank.getName().replace('.gz', '')
    template(task.ext.template)
}


/*
process primer_query {
    /*
    Query a set of PCR primers (FASTA) against annotated assembly using BLAST
    /

    cpus cpus
    tag "${sample} - ${reference}"
    publishDir "${outdir}/${sample}/primers", mode: 'copy', overwrite: true

    input:
    set val(sample), file(blast_db) from PRIMER_QUERY
    each primer from PRIMERS

    output:
    file("${sample}*.txt"

    shell:
    template(task.ext.template)
}
*/

workflow.onComplete {
    if (workflow.success == true && params.keep_cache == false) {
        clean_cache()
    }
    println """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    Launch Dir  : ${workflow.launchDir}
    Working Dir : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Command line: ${workflow.commandLine}
    Resumed?    : ${workflow.resume}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

// Utility functions
def print_usage(full_usage) {
    usage_text = ""
    if (params.full_usage) {
        // Print Full Usage
        usage_text = new File("$baseDir/configs/usage_full.txt").text
    } else {
        // Print basic usage
        usage_text = new File("$baseDir/configs/usage_basic.txt").text
    }
    log.info"""
        ${PROGRAM_NAME} v${VERSION}
        ${usage_text}
    """.stripIndent()
    clean_cache()
    exit 0
}


def print_version() {
    println(PROGRAM_NAME + ' ' + VERSION)
    clean_cache()
    exit 0
}

def print_example_fastqs() {
    log.info 'Printing example input for "--fastqs"'
    log.info ''
    log.info 'sample\tr1\tr2'
    log.info 'test001\t/path/to/fastqs/test_R1.fastq.gz\t/path/to/fastqs/test_R2.fastq.gz'
    log.info 'test002\t/path/to/fastqs/test.fastq.gz\t'
    clean_cache()
    exit 0
}

def print_check_fastqs(fastq_input) {
    log.info 'Printing what would have been processed. Each line consists of an array of'
    log.info 'three elements: [SAMPLE_NAME, IS_SINGLE_END, [FASTQ_1, FASTQ_2]]'
    log.info ''
    log.info 'Found:'
    create_fastq_channel(fastq_input).println()
    clean_cache()
    exit 0
}

def print_available_databases(database) {
    exit_code = 0
    if (database) {
        if (file("${database}/summary.json").exists()) {
            available_databases = read_database_summary(database)
            log.info 'Printing the available pre-configured databases.'
            log.info "Database Location (--database): ${database}"
            log.info ''
            if (available_databases.size() > 0) {
                available_databases.each { key, value ->
                    if (key == 'ariba') {
                        log.info "${key.capitalize()}"
                        value.each {
                            log.info "\tFound ${it}"
                        }
                    } else if (key == 'minmer') {
                        log.info "Minmer Sketches"
                        value.each {
                            log.info "\tFound ${it})"
                        }
                    } else {
                        log.info "${key.capitalize().replace('-', ' ')} (use --organism \"${key}\")"
                        value.each {
                            log.info "\tFound ${it}"
                        }
                    }
                    log.info ''
                }
            }
        } else {
            log.info "Please verify the PATH is correct and ${database}/summary.json" +
                     " exists, if not try rerunning 'setup-databases.py'."
            exit_code = 1
        }
    } else {
        log.info "Please use '--database' to specify the path to pre-built databases."
        exit_code = 1
    }
    clean_cache()
    exit exit_code
}

def read_database_summary(database) {
    slurp = new JsonSlurper()
    return slurp.parseText(file("${database}/summary.json").text)
}

def check_organism_databases(database, organism) {
    /* Check for available organism specific databases */
    organism_databases = []
    files = [
        // MLST
        "${database}/${organism}/mlst/ariba/ref_db/00.auto_metadata.tsv",
        "${database}/${organism}/mlst/blast/profile.txt",

        // Prokka
        "${database}/${organism}/prokka/proteins.faa"
    ]
    files.each {
        if (file(it).exists()) {
            organism_databases << it
        }
    }
    return organism_databases
}


def check_ariba_databases(database) {
    /* Check for available Ariba databases */
    ariba_directory = new File("${database}/ariba/")
    ariba_databases = []
    ariba_directory.eachFile(FileType.DIRECTORIES) {
        ariba_databases << it.name
    }
    return ariba_databases
}


def check_minmers(database) {
    /* Check for available minmer sketches */
    minmers = []
    sketches = [
        // Mash
        'refseq-k21-s1000.msh',
        // Sourmash
        'genbank-k21.json.gz', 'genbank-k31.json.gz', 'genbank-k51.json.gz'
    ]
    sketches.each {
        if (file("${database}/minmer/${it}").exists()) {
            minmers << "${database}/minmer/${it}"
        }
    }

    return minmers
}

def clean_cache() {
    // No need to resume completed run so remove cache.
    if (params.keep_cache == false) {
        file('./work/').deleteDir()
        file('./.nextflow/').deleteDir()
        file('.nextflow.log').delete()
    }
}

def is_positive_integer(value, name) {
    is_positive = true
    if (value.getClass() == Integer) {
        if (value < 0) {
            log.info('Invalid input (--'+ name +'), "' + value + '"" is not a positive integer.')
            is_positive = false
        }
    } else {
        if (!value.isInteger()) {
            log.info('Invalid input (--'+ name +'), "' + value + '"" is not numeric.')
            is_positive = false
        } else if (value.toInteger() < 0) {
            log.info('Invalid input (--'+ name +'), "' + value + '"" is not a positive integer.')
            is_positive = false
        }
    }
    return is_positive
}

def check_input_params() {
    error = false
    missing_requirement = false
    if (!params.fastqs) {
        log.info('Missing required "--fastqs" input')
        error = true
    } else if (!file(params.fastqs).exists()) {
        log.info('Invalid input (--fastqs), please verify "' + params.fastqs + '"" exists.')
        error = true
    }

    if (!is_positive_integer(params.max_cpus, 'max_cpus')) {
        error = true
    }

    if (!is_positive_integer(params.cpus, 'cpus')) {
        error = true
    }

    if (params.genome_size) {
        if (!['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
            if (!is_positive_integer(params.genome_size, 'genome_size')) {
                error = true
            }
        }
    }

    if (params.adapters) {
        if (!file(params.adapters).exists()) {
            log.info('Invalid input (--adapters), please verify "' + params.adapters + '"" exists.')
            error = true
        }
    }
    if (params.phix) {
        if (!file(params.phix).exists()) {
            log.info('Invalid input (--phix), please verify "' + params.phix + '"" exists.')
            error = true
        }
    }

    if (error) {
        log.info('See --help for more information')
        exit 1
    }
}

def process_tsv(line) {
    /* Parse line and determine if single end or paired reads*/
    if (line.r2) {
        // Paired
        return tuple(line.sample, false, [file(line.r1), file(line.r2)])
    } else {
        // Single End
        return tuple(line.sample, true, [file(line.r1)])
    }
}

def create_fastq_channel(fastq_input) {
    return Channel.fromPath( file(fastq_input) )
            .splitCsv(header: true, sep: '\t')
            .map { row -> process_tsv(row) }
}

def check_input_fastqs(fastq_input) {
    /* Read through --fastqs and verify each input exists. */
    samples = [:]
    error = false
    has_valid_header = false
    line = 1
    file(fastq_input).splitEachLine('\t') { cols ->
        if (line == 1) {
            if (cols[0] == 'sample' && cols[1] == 'r1' && cols[2] == 'r2') {
                has_valid_header = true
            }
        } else {
            if (samples.containsKey(cols[0])) {
                samples[cols[0]] = samples[cols[0]] + 1
            } else {
                samples[cols[0]] = 1
            }
            if (cols[1]) {
                if (!file(cols[1]).exists()) {
                    log.info "LINE " + line + ':ERROR: Please verify ' + cols[1]+ ' exists, and try again'
                    error = true
                }
            }
            if (cols[2]) {
                if (!file(cols[2]).exists()) {
                    log.info "LINE " + line + ':ERROR: Please verify ' + cols[2]+ ' exists, and try again'
                    error = true
                }
            }
        }

        line = line + 1
    }

    samples.each{ sample, count ->
        if (count > 1) {
            error = true
            log.info 'Sample name "'+ sample +'" is not unique, please revise sample names'
        }
    }

    if (!has_valid_header) {
        error = true
        log.info 'The header line (line 1) does not follow expected structure.'
    }

    if (error) {
        log.info 'Verify sample names are unique and/or FASTQ paths are correct'
        log.info 'See "--example_fastqs" for an example'
        log.info 'Exiting'
        exit 1
    }
}

def get_absolute_path(file_path) {
    // Thanks Fabian Steeg
    // https://stackoverflow.com/questions/3204955/converting-relative-paths-to-absolute-paths
    File file_obj = new File("${workflow.launchDir}/${file_path}");
    String absolute_path = file_obj.getCanonicalPath(); // may throw IOException

    return absolute_path
}

def print_database_info(database_list, database_info) {
    log.info "Found ${database_list.size()} ${database_info}"
    database_list.each {
        log.info "\t${it}"
    }
}
