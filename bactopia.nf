#! /usr/bin/env nextflow
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine
import java.nio.file.Path
import java.nio.file.Paths
PROGRAM_NAME = 'bactopia'
VERSION = '0.0.1'
if (params.help || params.help_all) print_usage();
if (workflow.commandLine.endsWith(workflow.scriptName)) print_usage();
if (params.example_fastqs) print_example_fastqs();
if (params.version) print_version();
check_input_params()
check_input_fastqs(params.fastqs)
if (params.check_fastqs) print_check_fastqs(params.fastqs);
if (params.available_datasets) print_available_datasets(params.dataset)

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
BLAST_FASTAS = []
MAPPING_FASTAS = []
PLASMID_BLASTDB = []
PROKKA_PROTEINS = null
species_genome_size = ['min': 0, 'median': 0, 'mean': 0, 'max': 0]
if (params.dataset) {
    dataset_path = get_canonical_path(params.dataset)
    available_datasets = read_dataset_summary(dataset_path)

    available_datasets['ariba'].each {
        ARIBA_DATABASES << file("${dataset_path}/ariba/${it.name}")
    }
    print_dataset_info(ARIBA_DATABASES, "ARIBA datasets")

    available_datasets['minmer']['sketches'].each {
        MINMER_DATABASES << file("${dataset_path}/minmer/${it}")
    }
    MINMER_DATABASES << file("${dataset_path}/plasmid/${available_datasets['plasmid']['sketches']}")
    print_dataset_info(MINMER_DATABASES, "minmer sketches/signatures")

    PLASMID_BLASTDB = tuple(file("${dataset_path}/plasmid/${available_datasets['plasmid']['blastdb']}*"))
    print_dataset_info(PLASMID_BLASTDB, "PLSDB (plasmid) BLAST files")

    if (params.species) {
        if (available_datasets['species-specific'].containsKey(params.species)) {
            species_db = available_datasets['species-specific'][params.species]
            species_genome_size = species_db['genome_size']

            prokka = "${dataset_path}/${species_db['annotation']['proteins']}"
            if (file(prokka).exists()) {
                PROKKA_PROTEINS = file(prokka)
                log.info "Found Prokka proteins file"
                log.info "\t${PROKKA_PROTEINS}"

            }
            species_db['mlst'].each { key, val ->
                if (key != "last_updated") {
                    if (file("${dataset_path}/${val}").exists()) {
                        MLST_DATABASES << file("${dataset_path}/${val}")
                    }
                }
            }
            print_dataset_info(MLST_DATABASES, "MLST datasets")

            file("${dataset_path}/${species_db['optional']['reference-genomes']}").list().each() {
                REFERENCES << file("${dataset_path}/${species_db['optional']['reference-genomes']}/${it}")
            }
            print_dataset_info(REFERENCES, "reference genomes")

            file("${dataset_path}/${species_db['optional']['insertion-sequences']}").list().each() {
                INSERTIONS << file("${dataset_path}/${species_db['optional']['insertion-sequences']}/${it}")
            }
            print_dataset_info(INSERTIONS, "insertion sequence FASTAs")

            file("${dataset_path}/${species_db['optional']['mapping-sequences']}").list().each() {
                MAPPING_FASTAS << file("${dataset_path}/${species_db['optional']['mapping-sequences']}/${it}")
            }
            print_dataset_info(MAPPING_FASTAS, "FASTAs to align reads against")

            // BLAST Related
            species_db['optional']['blast'].each() {
                temp_path = "${dataset_path}/${it}"
                file(temp_path).list().each() {
                    BLAST_FASTAS << file("${temp_path}/${it}")
                }
            }
            print_dataset_info(BLAST_FASTAS, "FASTAs to query with BLAST")
        } else {
            log.info "Species '${params.species}' not available, please check spelling or use '--available_datasets' " +
                     "to verify the dataset has been set up. Exiting"
            exit 1
        }
    } else {
        log.info "--species not given, skipping the following processes (analyses):"
        log.info "\tsequence_type"
        log.info "\tcall_variants"
        log.info "\tinsertion_sequence_query"
        log.info "\tprimer_query"
        if (['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
            log.info "Asked for genome size '${params.genome_size}' which requires a " +
                     "species to be given. Please give a species or specify " +
                     "a valid genome size. Exiting"
            exit 1
        }
    }
} else {
    log.info "--dataset not given, skipping the following processes (analyses):"
    log.info "\tsequence_type"
    log.info "\tariba_analysis"
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
                                                         INSERTION_SEQUENCES, CALL_VARIANTS, MAPPING_QUERY

    shell:
    fq2 = single_end == true ? "" : fq[1]
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
    set val(sample), file("${sample}.fna.gz") into ANNOTATION, MAKE_BLASTDB

    shell:
    opts = params.shovill_opts ? "--opts '${params.shovill_opts}'" : ""
    kmers = params.shovill_kmers ? "--kmers '${params.shovill_kmers}'" : ""
    nostitch = params.nostitch ? "--nostitch" : ""
    nocorr = params.nocorr ? "--nocorr" : ""
    template(task.ext.template)
}


process make_blastdb {
    /* Create a BLAST database of the assembly using BLAST */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), file(fasta) from MAKE_BLASTDB

    output:
    set val(sample), file("blastdb/*") into BLAST_QUERY

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
    set val(sample), file("annotation/*.ffn.gz") into PLASMID_BLAST

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    proteins = PROKKA_PROTEINS ? "--proteins ${PROKKA_PROTEINS}" : ""
    genus = PROKKA_PROTEINS ? params.species.split('-')[0].capitalize() : "Genus"
    species = PROKKA_PROTEINS ? params.species.split('-')[1] : "species"
    addgenes = params.addgenes ? "--addgenes" : ""
    addmrna = params.addmrna ? "--addmrna" : ""
    rawproduct = params.rawproduct ? "--rawproduct" : ""
    cdsrnaolap = params.cdsrnaolap ? "--cdsrnaolap" : ""
    norrna = params.norrna ? "--norrna" : ""
    notrna = params.notrna ? "--notrna" : ""
    rnammer = params.rnammer ? "--rnammer" : ""
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
    cpus { task.attempt > 1 ? 1 : cpus }
    errorStrategy 'retry'
    maxRetries 5
    tag "${sample} - ${method}"
    publishDir "${outdir}/${sample}/mlst", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from SEQUENCE_TYPE
    file(assembly) from SEQUENCE_TYPE_ASSEMBLY
    each dataset from MLST_DATABASES

    when:
    dataset.contains('blast') || (dataset.contains('ariba') && single_end == false)

    output:
    file "${method}/*"

    shell:
    method = dataset.contains('blast') ? 'blast' : 'ariba'
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    template(task.ext.template)
}


process ariba_analysis {
    /* Run reads against all available (if any) ARIBA datasets */
    cpus { task.attempt > 1 ? 1 : cpus }
    errorStrategy 'retry'
    maxRetries 5
    tag "${sample} - ${dataset_name}"
    publishDir "${outdir}/${sample}/ariba", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from ARIBA_ANALYSIS
    each dataset from ARIBA_DATABASES

    output:
    file "${dataset_name}/*"

    when:
    single_end == false

    shell:
    dataset_name = file(dataset).getName()
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
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
    tag "${sample} - ${dataset_name}"
    publishDir "${outdir}/${sample}/minmers", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from MINMER_QUERY
    file(sourmash) from QUERY_SOURMASH
    each file(dataset) from MINMER_DATABASES

    output:
    file("${sample}*.txt")

    shell:
    dataset_name = dataset.getName()
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
    each file(reference) from REFERENCES

    output:
    file("${reference_name}/*")

    shell:
    reference_name = reference.getSimpleName()
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
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
    each file(insertion_fasta) from INSERTIONS

    output:
    file("insertion-sequences/*")

    when:
    single_end == false

    shell:
    insertion_name = insertion_fasta.getSimpleName()
    gunzip_genbank = genbank.getName().replace('.gz', '')
    all = params.ismap_all ? "--a" : ""
    template(task.ext.template)
}


process plasmid_blast {
    /*
    BLAST a set of predicted genes against the PLSDB BALST database.
    */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/plasmid", mode: 'copy', overwrite: true

    input:
    set val(sample), file(genes) from PLASMID_BLAST
    file(blastdb_files) from Channel.from(PLASMID_BLASTDB).toList()

    output:
    file("${sample}-plsdb.txt")

    shell:
    blastdb = blastdb_files[0].getBaseName()
    template(task.ext.template)
}


process blast_query {
    /*
    Query a FASTA files against annotated assembly using BLAST
    */

    cpus cpus
    tag "${sample} - ${query.getName()}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), file(blastdb) from BLAST_QUERY
    each file(query) from BLAST_FASTAS

    output:
    file("blast/*")

    shell:
    query_name = query.getSimpleName()
    template(task.ext.template)
}


process mapping_query {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */

    cpus cpus
    tag "${sample} - ${query.getName()}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from MAPPING_QUERY
    each file(query) from MAPPING_FASTAS

    output:
    file("mapping/*")

    shell:
    query_name = query.getSimpleName()
    bwa_mem_opts = params.bwa_mem_opts ? params.bwa_mem_opts : ""
    bwa_aln_opts = params.bwa_aln_opts ? params.bwa_aln_opts : ""
    bwa_samse_opts = params.bwa_samse_opts ? params.bwa_samse_opts : ""
    bwa_sampe_opts = params.bwa_sampe_opts ? params.bwa_sampe_opts : ""
    template(task.ext.template)
}


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
def print_usage() {
    usage_text = ""
    if (params.help_all) {
        // Print Full Usage
        usage_text = basic_help()
    } else {
        // Print basic usage
        usage_text = basic_help()
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

def print_available_datasets(dataset) {
    exit_code = 0
    if (dataset) {
        if (file("${dataset}/summary.json").exists()) {
            available_datasets = read_dataset_summary(dataset)
            log.info 'Printing the available pre-configured dataset.'
            log.info "Database Location (--dataset): ${dataset}"
            log.info ''
            if (available_datasets.size() > 0) {
                available_datasets.each { key, value ->
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
                        log.info "${key.capitalize().replace('-', ' ')} (use --species \"${key}\")"
                        value.each {
                            log.info "\tFound ${it}"
                        }
                    }
                    log.info ''
                }
            }
        } else {
            log.info "Please verify the PATH is correct and ${dataset}/summary.json" +
                     " exists, if not try rerunning 'setup-public-datasets.py'."
            exit_code = 1
        }
    } else {
        log.info "Please use '--dataset' to specify the path to pre-built datasets."
        exit_code = 1
    }
    clean_cache()
    exit exit_code
}

def read_dataset_summary(dataset) {
    slurp = new JsonSlurper()
    return slurp.parseText(file("${dataset}/summary.json").text)
}

def check_species_datasets(dataset, species) {
    /* Check for available species specific datasets */
    species_datasets = []
    files = [
        // MLST
        "${dataset}/${species}/mlst/ariba/ref_db/00.auto_metadata.tsv",
        "${dataset}/${species}/mlst/blast/profile.txt",

        // Prokka
        "${dataset}/${species}/prokka/proteins.faa"
    ]
    files.each {
        if (file(it).exists()) {
            species_datasets << it
        }
    }
    return species_datasets
}


def check_ariba_datasets(dataset) {
    /* Check for available Ariba datasets */
    ariba_directory = new File("${dataset}/ariba/")
    ariba_datasets = []
    ariba_directory.eachFile(FileType.DIRECTORIES) {
        ariba_datasets << it.name
    }
    return ariba_datasets
}


def check_minmers(dataset) {
    /* Check for available minmer sketches */
    minmers = []
    sketches = [
        // Mash
        'refseq-k21-s1000.msh',
        // Sourmash
        'genbank-k21.json.gz', 'genbank-k31.json.gz', 'genbank-k51.json.gz'
    ]
    sketches.each {
        if (file("${dataset}/minmer/${it}").exists()) {
            minmers << "${dataset}/minmer/${it}"
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


def process_csv(line) {
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
        .map { row -> process_csv(row) }
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
    String absolute_path = file_obj.getAbsolutePath(); // may throw IOException

    return absolute_path
}

def get_canonical_path(file_path) {
    // Thanks Fabian Steeg
    // https://stackoverflow.com/questions/3204955/converting-relative-paths-to-absolute-paths
    File file_obj = new File("${workflow.launchDir}/${file_path}");
    String canonical_path = file_obj.getCanonicalPath(); // may throw IOException

    return canonical_path
}

def get_real_path(file_path) {
    Path real_path = Paths.get(file_path)
    println real_path.toRealPath()
    return real_path.toRealPath()
}


def print_dataset_info(dataset_list, dataset_info) {
    log.info "Found ${dataset_list.size()} ${dataset_info}"
    dataset_list.each {
        log.info "\t${it}"
    }
}

def basic_help() {
    genome_size = params.genome_size ? params.genome_size : "Mash Estimate"
    return """
    Required Parameters:
        --fastqs STR            An input file containing the sample name and
                                    absolute paths to FASTQs to process

    Public Dataset Parameters:
        --datasets DIR          The path to available public datasets that have
                                    already been set up
        --species STR           Determines which species-specific dataset to use
                                    for the input sequencing

    Optional Parameters:
        --coverage INT          Reduce samples to a given coverage
                                    (Default: ${params.coverage}x)

        --genome_size INT       Expected genome size (bp) for all samples
                                    (Default: ${genome_size})

        --outdir DIR            Directory to write results to (Default ${params.outdir})

        --max_cpus INT          The maximum number of processors this workflow
                                    should have access to at any given moment
                                    (Default: ${params.max_cpus})

        --cpus INT              Number of processors made available to a single
                                    process. If greater than "--max_cpus" it will
                                    be set equal to "--max_cpus"
                                    (Default: ${params.cpus})

    Useful Parameters:
        --available_datasets    Print a list of available datasets found based
                                    on location given by "--datasets"
        --example_fastqs        Print example of expected input for FASTQs file

        --check_fastqs          Verify "--fastqs" produces the expected inputs

        --keep_cache            Keeps 'work' and '.nextflow' logs, default is
                                    to delete on successful completion

        --version               Print workflow version information
        --help                  Show this message and exit
        --help_all              Show a complete list of adjustable parameters
    """

}
