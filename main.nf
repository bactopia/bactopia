#! /usr/bin/env nextflow
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine
import groovy.util.FileNameByRegexFinder
import java.nio.file.Path
import java.nio.file.Paths
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version

// Validate parameters
if (params.help || params.help_all || params.conda_help) print_usage();
if (params.nfdir) print_basedir();
if (workflow.commandLine.trim().endsWith(workflow.scriptName)) print_usage();
if (params.example_fastqs) print_example_fastqs();
if (params.version) print_version();
fastq_type = check_input_params()
check_input_fastqs(params.fastqs, fastq_type)
if (params.check_fastqs) print_check_fastqs(params.fastqs, fastq_type);
if (params.available_datasets) print_available_datasets(params.dataset)

// Setup output directories
outdir = params.outdir ? params.outdir : './'

// Setup some defaults
log.info "${PROGRAM_NAME} - ${VERSION}"
ARIBA_DATABASES = []
MINMER_DATABASES = []
MLST_DATABASES = []
REFERENCES = []
INSERTIONS = []
PRIMERS = []
BLAST_FASTAS = []
MAPPING_FASTAS = []
PLASMID_BLASTDB = []
PROKKA_PROTEINS = file('EMPTY')
REFSEQ_SKETCH = []
REFSEQ_SKETCH_FOUND = false
SPECIES = format_species(params.species)
setup_datasets()


process gather_fastqs {
    /* Gather up input FASTQs for analysis. */

    input:
    set val(sample), val(single_end), file(fq) from create_fastq_channel(params.fastqs, fastq_type)

    output:
    set val(sample), val(single_end),
        file("fastqs/${sample}*.fastq.gz") into FASTQ_PE_STATUS

    shell:
    template(task.ext.template)
}

process fastq_status {
    /* Determine if FASTQs are PE or SE, and if they meet minimum basepair/read counts. */
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: '*.txt'

    input:
    set val(sample), val(single_end), file(fq) from FASTQ_PE_STATUS

    output:
    file "*-error.txt" optional true
    set val(sample), val(single_end), file("fastqs/${sample}*.fastq.gz") optional true into ESTIMATE_GENOME_SIZE

    shell:
    single_end = fq[1] == null ? true : false
    template(task.ext.template)
}

process estimate_genome_size {
    /* Estimate the input genome size if not given. */
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: '*.txt'

    input:
    set val(sample), val(single_end), file(fq) from ESTIMATE_GENOME_SIZE

    output:
    file "${sample}-genome-size-error.txt" optional true
    file("${sample}-genome-size.txt") optional true
    set val(sample), val(single_end), file(fq), file("${sample}-genome-size.txt") optional true into QC_READS, QC_ORIGINAL_SUMMARY

    shell:
    template(task.ext.template)
}

process qc_reads {
    /* Cleanup the reads using Illumina-Cleanup */
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "quality-control/*"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "*error.txt"

    input:
    set val(sample), val(single_end), file(fq), file(genome_size) from QC_READS

    output:
    file "*-error.txt" optional true
    file "quality-control/*"
    set val(sample), val(single_end),
        file("quality-control/${sample}*.fastq.gz") optional true into COUNT_31MERS, ARIBA_ANALYSIS,
                                                                       MINMER_SKETCH, CALL_VARIANTS,
                                                                       MAPPING_QUERY
    set val(sample), val(single_end),
        file("quality-control/${sample}*.fastq.gz"),
        file(genome_size) optional true into ASSEMBLY

    set val(sample), val(single_end),
        file("quality-control/${sample}*.{fastq,error-fq}.gz"),
        file(genome_size) optional true into QC_FINAL_SUMMARY

    shell:
    adapters = params.adapters ? file(params.adapters) : 'adapters'
    phix = params.phix ? file(params.phix) : 'phix'
    template(task.ext.template)
}


process qc_original_summary {
    /* Run FASTQC on the input FASTQ files. */
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    set val(sample), val(single_end), file(fq), file(genome_size) from QC_ORIGINAL_SUMMARY

    output:
    file "quality-control/*"

    when:
    params.dry_run == false

    shell:
    template(task.ext.template)
}

process qc_final_summary {
    /* Run FASTQC on the cleaned up FASTQ files. */
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "quality-control/*"

    input:
    set val(sample), val(single_end), file(fq), file(genome_size) from QC_FINAL_SUMMARY

    output:
    file "quality-control/*"

    when:
    params.dry_run == false

    shell:
    template(task.ext.template)
}


process assemble_genome {
    /* Assemble the genome using Shovill, SKESA is used by default */
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "assembly/*"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "${sample}-assembly-error.txt"

    input:
    set val(sample), val(single_end), file(fq), file(genome_size) from ASSEMBLY

    output:
    file "assembly/*"
    file "${sample}-assembly-error.txt" optional true
    set val(sample), val(single_end), file(fq), file("assembly/${sample}.{fna,fna.gz}") optional true into ANNOTATION, MAKE_BLASTDB, SEQUENCE_TYPE

    shell:
    shovill_ram = task.memory.toString().split(' ')[0]
    opts = params.shovill_opts ? "--opts '${params.shovill_opts}'" : ""
    kmers = params.shovill_kmers ? "--kmers '${params.shovill_kmers}'" : ""
    nostitch = params.nostitch ? "--nostitch" : ""
    nocorr = params.nocorr ? "--nocorr" : ""
    template(task.ext.template)
}


process make_blastdb {
    /* Create a BLAST database of the assembly using BLAST */
    tag "${sample}"
    publishDir "${outdir}/${sample}/blast", mode: 'copy', overwrite: params.overwrite, pattern: "blastdb/*"

    input:
    set val(sample), val(single_end), file(fq), file(fasta) from MAKE_BLASTDB

    output:
    file("blastdb/*")
    set val(sample), file("blastdb/*") into BLAST_QUERY

    shell:
    template(task.ext.template)
}


process annotate_genome {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "annotation/${sample}*"

    input:
    set val(sample), val(single_end), file(fq), file(fasta) from ANNOTATION
    file prokka_proteins from PROKKA_PROTEINS

    output:
    file "annotation/${sample}*"
    set val(sample), file("annotation/${sample}.{ffn,ffn.gz}") optional true into PLASMID_BLAST
    set val(sample), val(single_end), file(fq), file("annotation/${sample}.{gbk,gbk.gz}") optional true into INSERTION_SEQUENCES
    set val(sample),
        file("annotation/${sample}.{ffn,ffn.gz}"),
        file("annotation/${sample}.{faa,faa.gz}") optional true into ANTIMICROBIAL_RESISTANCE

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    genus = "Genus"
    species = "species"
    proteins = ""
    if (prokka_proteins.getName() != 'EMPTY') {
        proteins = "--proteins ${prokka_proteins}"
        if (SPECIES.contains("-")) {
            genus = SPECIES.split('-')[0].capitalize()
            species = SPECIES.split('-')[1]
        } else {
            genus = SPECIES.capitalize()
            species = "spp."
        }
    }
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
    tag "${sample}"
    publishDir "${outdir}/${sample}/kmers", mode: 'copy', overwrite: params.overwrite, pattern: "*.ctx"

    input:
    set val(sample), val(single_end), file(fq) from COUNT_31MERS

    output:
    file "${sample}.ctx"

    shell:
    m = task.memory.toString().split(' ')[0].toInteger() * 1000
    template(task.ext.template)
}


process sequence_type {
    /* Determine MLST types using ARIBA and BLAST */
    tag "${sample} - ${method}"
    publishDir "${outdir}/${sample}/mlst", mode: 'copy', overwrite: params.overwrite, pattern: "${method}/*"

    input:
    set val(sample), val(single_end), file(fq), file(assembly) from SEQUENCE_TYPE
    each file(dataset) from MLST_DATABASES

    when:
    dataset =~ /.*blast.*/ || (dataset =~ /.*ariba.*/ && single_end == false)

    output:
    file "${method}/*"

    shell:
    method = dataset =~ /.*blastdb.*/ ? 'blast' : 'ariba'
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    noclean = params.ariba_no_clean ? "--noclean" : ""
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    template(task.ext.template)
}


process ariba_analysis {
    /* Run reads against all available (if any) ARIBA datasets */
    tag "${sample} - ${dataset_name}"
    publishDir "${outdir}/${sample}/ariba", mode: 'copy', overwrite: params.overwrite, pattern: "${dataset_name}/*"

    input:
    set val(sample), val(single_end), file(fq) from ARIBA_ANALYSIS
    each file(dataset) from ARIBA_DATABASES

    output:
    file "${dataset_name}/*"

    when:
    single_end == false || params.dry_run == true

    shell:
    dataset_tarball = file(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    noclean = params.ariba_no_clean ? "--noclean" : ""
    template(task.ext.template)
}


process minmer_sketch {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31) and
    Sourmash (k=21,31,51)
    */
    tag "${sample}"
    publishDir "${outdir}/${sample}/minmers", mode: 'copy', overwrite: params.overwrite, pattern: "*.{msh,sig}"

    input:
    set val(sample), val(single_end), file(fq) from MINMER_SKETCH

    output:
    file("${sample}*.{msh,sig}")
    set val(sample), val(single_end), file(fq), file("${sample}.sig") into MINMER_QUERY
    set val(sample), val(single_end), file(fq), file("${sample}-k31.msh") into DOWNLOAD_REFERENCES

    shell:
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    template(task.ext.template)
}


process minmer_query {
    /*
    Query minmer sketches against pre-computed RefSeq (Mash, k=21) and
    GenBank (Sourmash, k=21,31,51)
    */
    tag "${sample} - ${dataset_name}"
    publishDir "${outdir}/${sample}/minmers", mode: 'copy', overwrite: params.overwrite, pattern: "*.txt"

    input:
    set val(sample), val(single_end), file(fq), file(sourmash) from MINMER_QUERY
    each file(dataset) from MINMER_DATABASES

    output:
    file "*.txt"

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
    tag "${sample} - ${reference_name}"
    publishDir "${outdir}/${sample}/variants/user", mode: 'copy', overwrite: params.overwrite, pattern: "${reference_name}/*"

    input:
    set val(sample), val(single_end), file(fq) from CALL_VARIANTS
    each file(reference) from REFERENCES

    output:
    file("${reference_name}/*")

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getSimpleName()
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    template(task.ext.template)
}


process download_references {
    /*
    Download the nearest RefSeq genomes (based on Mash) to have variants called against.

    Exitcode 75 is due to being unable to download from NCBI (e.g. FTP down at the time)
    Downloads will be attempted 300 times total before giving up. On failure to download
    variants will not be called against the nearest completed genome.
    */
    tag "${sample} - ${params.max_references} reference(s)"
    publishDir "${outdir}/${sample}/variants/auto", mode: 'copy', overwrite: params.overwrite, pattern: 'mash-dist.txt'

    input:
    set val(sample), val(single_end), file(fq), file(sample_sketch) from DOWNLOAD_REFERENCES
    file(refseq_sketch) from REFSEQ_SKETCH

    output:
    set val(sample), val(single_end), file(fq), file("genbank/*.gbk") optional true into CALL_VARIANTS_AUTO
    file("mash-dist.txt")

    when:
    REFSEQ_SKETCH_FOUND == true || params.dry_run == true

    shell:
    tie_break = params.random_tie_break ? "--random_tie_break" : ""
    total = params.max_references
    template(task.ext.template)
}


process call_variants_auto {
    /*
    Identify variants (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample} - ${reference_name}"
    publishDir "${outdir}/${sample}/variants/auto", mode: 'copy', overwrite: params.overwrite, pattern: "${reference_name}/*"

    input:
    set val(sample), val(single_end), file(fq), file(reference) from create_reference_channel(CALL_VARIANTS_AUTO)

    output:
    file("${reference_name}/*")

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getBaseName().split('-')[1].split(/\./)[0]
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    template(task.ext.template)
}


process update_antimicrobial_resistance {
    /*
    Update amrfinders database a single time.
    */

    input:
    val x from Channel.from(1)

    output:
    val(updated) into AMR_UPDATED

    shell:
    updated = true
    template(task.ext.template)
}


process antimicrobial_resistance {
    /*
    Query nucleotides and proteins (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "${amrdir}/*"

    input:
    set val(sample), file(genes), file(proteins) from ANTIMICROBIAL_RESISTANCE

    output:
    file("${amrdir}/*")

    shell:
    amrdir = "antimicrobial_resistance"
    plus = params.amr_plus ? "--plus" : ""
    report_common = params.amr_report_common ? "--report_common" : ""
    organism_gene = ""
    organism_protein = ""
    if (params.amr_organism) {
        organism_gene = "-O ${params.amr_organism} --point_mut_all ${amrdir}/${sample}-gene-point-mutations.txt"
        organism_protein = "-O ${params.amr_organism} --point_mut_all ${amrdir}/${sample}-protein-point-mutations.txt"
    }
    template(task.ext.template)
}


process insertion_sequences {
    /*
    Query a set of insertion sequences (FASTA) against annotated GenBank file
    using ISMapper.
    */
    tag "${sample} - ${insertion_name}"
    publishDir "${outdir}/${sample}/", mode: 'copy', overwrite: params.overwrite, pattern: "insertion-sequences/*"

    input:
    set val(sample), val(single_end), file(fq), file(genbank) from INSERTION_SEQUENCES
    each file(insertion_fasta) from INSERTIONS

    output:
    file("insertion-sequences/*")

    when:
    single_end == false || params.dry_run == true

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
    tag "${sample}"
    publishDir "${outdir}/${sample}/blast", mode: 'copy', overwrite: params.overwrite, pattern: "*.{txt,txt.gz}"

    input:
    set val(sample), file(genes) from PLASMID_BLAST
    file(blastdb_files) from Channel.from(PLASMID_BLASTDB).toList()

    output:
    file("${sample}-plsdb.{txt,txt.gz}")

    when:
    PLASMID_BLASTDB.isEmpty() == false

    shell:
    blastdb = blastdb_files[0].getBaseName()
    template(task.ext.template)
}


process blast_query {
    /*
    Query a FASTA files against annotated assembly using BLAST
    */
    tag "${sample} - ${query.getName()}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "blast/*"

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
    tag "${sample} - ${query.getName()}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: params.overwrite, pattern: "mapping/*"

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
    f_value = params.keep_unmapped_reads ? '-F 0' : '-F 4'
    template(task.ext.template)
}


workflow.onComplete {
    workDir = new File("${workflow.workDir}")
    workDirSize = toHumanString(workDir.directorySize())

    println """

    Bactopia Execution Summary
    ---------------------------
    Command Line    : ${workflow.commandLine}
    Resumed         : ${workflow.resume}

    Completed At    : ${workflow.complete}
    Duration        : ${workflow.duration}
    Success         : ${workflow.success}
    Exit Code       : ${workflow.exitStatus}
    Error Report    : ${workflow.errorReport ?: '-'}

    Launch Dir      : ${workflow.launchDir}
    Working Dir     : ${workflow.workDir} (Total Size: ${workDirSize})
    Working Dir Size: ${workDirSize}

    """
}


// Utility functions
def toHumanString(bytes) {
    // Thanks Niklaus
    // https://gist.github.com/nikbucher/9687112
    base = 1024L
    decimals = 3
    prefix = ['', 'K', 'M', 'G', 'T']
    int i = Math.log(bytes)/Math.log(base) as Integer
    i = (i >= prefix.size() ? prefix.size()-1 : i)
    return Math.round((bytes / base**i) * 10**decimals) / 10**decimals + prefix[i]
}

def print_version() {
    println(PROGRAM_NAME + ' ' + VERSION)
    exit 0
}


def print_basedir() {
    println("BACTOPIA_DIR = ${baseDir}")
    exit 0
}

def print_example_fastqs() {
    log.info 'Printing example input for "--fastqs"'
    log.info ''
    log.info 'sample\tr1\tr2'
    log.info 'test001\t/path/to/fastqs/test_R1.fastq.gz\t/path/to/fastqs/test_R2.fastq.gz'
    log.info 'test002\t/path/to/fastqs/test.fastq.gz\t'
    exit 0
}


def print_check_fastqs(fastq_input, fastq_type) {
    if (fastq_type == "fastqs") {
        log.info 'Printing what would have been processed. Each line consists of an array of'
        log.info 'three elements: [SAMPLE_NAME, IS_SINGLE_END, [FASTQ_1, FASTQ_2]]'
        log.info ''
        log.info 'Found:'
        create_fastq_channel(fastq_input).println()
        exit 0
    } else {
        log.error '"--check_fastqs" requires "--fastqs" to be set.'
        exit 1
    }
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
            log.error "Please verify the PATH is correct and ${dataset}/summary.json" +
                     " exists, if not try rerunning 'setup-datasets.py'."
            exit_code = 1
        }
    } else {
        log.error "Please use '--dataset' to specify the path to pre-built datasets."
        exit_code = 1
    }
    exit exit_code
}


def read_dataset_summary(dataset) {
    slurp = new JsonSlurper()
    return slurp.parseText(file("${dataset}/summary.json").text)
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


def setup_datasets() {
    species_genome_size = ['min': 0, 'median': 0, 'mean': 0, 'max': 0]
    if (params.datasets) {
        dataset_path = params.datasets
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

        if (SPECIES) {
            if (available_datasets.containsKey('species-specific')) {
                if (available_datasets['species-specific'].containsKey(SPECIES)) {
                    species_db = available_datasets['species-specific'][SPECIES]
                    species_genome_size = species_db['genome_size']

                    prokka = "${dataset_path}/${species_db['annotation']['proteins']}"
                    if (file(prokka).exists()) {
                        PROKKA_PROTEINS = file(prokka)
                        log.info "Found Prokka proteins file"
                        log.info "\t${PROKKA_PROTEINS}"
                    }

                    refseq_minmer = "${dataset_path}/${species_db['minmer']['mash']}"
                    if (file(refseq_minmer).exists()) {
                        REFSEQ_SKETCH = file(refseq_minmer)
                        REFSEQ_SKETCH_FOUND = true
                        log.info "Found Mash Sketch of RefSeq genomes"
                        log.info "\t${REFSEQ_SKETCH}"
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
            log.info "\tsequence_type"
            log.info "\tcall_variants"
            log.info "\tinsertion_sequence_query"
            log.info "\tprimer_query"
            if (['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
                log.error "Asked for genome size '${params.genome_size}' which requires a " +
                          "species to be given. Please give a species or specify " +
                          "a valid genome size. Exiting"
                exit 1
            }
        }
    } else if (params.dry_run) {
        log.info "--dry_run found, creating dummy data to be 'processed'."
        ARIBA_DATABASES << file('EMPTY.tar.gz')
        INSERTIONS << file('EMPTY.fasta')
        MAPPING_FASTAS << file('EMPTY.fasta')
        MLST_DATABASES << file('EMPTY_DB')
        REFSEQ_SKETCH = file('EMPTY.msh')
    } else {
        log.info "--datasets not given, skipping the following processes (analyses):"
        log.info "\tsequence_type"
        log.info "\tariba_analysis"
        log.info "\tminmer_query"
        log.info "\tplasmid_blast"
        log.info "\tcall_variants"
        log.info "\tinsertion_sequence_query"
        log.info "\tprimer_query"
    }

    if (params.disable_auto_variants) {
        REFSEQ_SKETCH_FOUND = false
    }


    log.info "\nIf something looks wrong, now's your chance to back out (CTRL+C 3 times). "
    log.info "Sleeping for ${params.sleep_time} seconds..."
    sleep(params.sleep_time * 1000)
}


def is_positive_integer(value, name) {
    error = 0
    if (value.getClass() == Integer) {
        if (value < 0) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not a positive integer.')
            error = 1
        }
    } else {
        if (!value.isInteger()) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not numeric.')
            error = 1
        } else if (value.toInteger() < 0) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not a positive integer.')
            error = 1
        }
    }
    return error
}


def file_exists(file_name, parameter) {
    if (!file(file_name).exists()) {
        log.error('Invalid input ('+ parameter +'), please verify "' + file_name + '" exists.')
        return 1
    }
    return 0
}


def check_unknown_params() {
    valid_params = []
    error = 0
    new File("${baseDir}/conf/params.config").eachLine { line ->
        if (line.contains("=")) {
            valid_params << line.trim().split(" ")[0]
        }
    }

    params.each { k,v ->
        if (!valid_params.contains(k)) {
            if (k != "container-path") {
                log.error("'--${k}' is not a known parameter")
                error = 1
            }
        }
    }

    return error
}


def check_input_params() {
    error = 0
    fastq_type = null

    // Check for unexpected paramaters
    error += check_unknown_params()

    if (params.fastqs) {
        error += file_exists(params.fastqs, '--fastqs')
        fastq_type = "fastqs"
    } else if  (params.R1 && params.R2 && params.sample) {
        error += file_exists(params.R1, '--R1')
        error += file_exists(params.R2, '--R2')
        fastq_type = "paired"
    } else if (params.SE && params.sample) {
        error += file_exists(params.SE, '--SE')
        fastq_type = "single"
    } else if (params.accessions) {
        error += file_exists(params.accessions, '--accessions')
        fastq_type = "ena_accessions"
    }else if (params.accession || params.dry_run) {
        fastq_type = "ena_accession"
    } else {
        log.error """
        One or more required parameters are missing, please check and try again.

        Required Parameters:
            ### For Procesessing Multiple Samples
            --fastqs STR            An input file containing the sample name and
                                        absolute paths to FASTQs to process

            ### For Processing A Single Sample
            --R1 STR                First set of reads for paired end in compressed (gzip)
                                        FASTQ format

            --R2 STR                Second set of reads for paired end in compressed (gzip)
                                        FASTQ format

            --SE STR                Single end set of reads in compressed (gzip) FASTQ format

            --sample STR            The name of the input sequences
            ### For Downloading from ENA
            --accessions            An input file containing ENA/SRA experiement accessions to
                                        be processed

            --accession             A single ENA/SRA Experiment accession to be processed
        """.stripIndent()
        error += 1
    }

    error += is_positive_integer(params.cpus, 'cpus')
    error += is_positive_integer(params.max_time, 'max_time')
    error += is_positive_integer(params.max_memory, 'max_memory')
    error += is_positive_integer(params.shovill_ram, 'shovill_ram')
    error += is_positive_integer(params.snippy_ram, 'snippy_ram')
    error += is_positive_integer(params.cortex_ram, 'cortex_ram')
    error += is_positive_integer(params.qc_ram, 'qc_ram')
    error += is_positive_integer(params.minmer_ram, 'minmer_ram')

    if (params.genome_size) {
        if (!['min', 'median', 'mean', 'max'].contains(params.genome_size)) {
            error += is_positive_integer(params.genome_size, 'genome_size')
        }
    }

    if (params.ftp_only && !params.use_ena) {
        log.error "'--ftp_only' must be used along side '--use_ena'"
        error += 1
    }

    if (params.adapters) {
        error += file_exists(params.adapters, '--adapters')
    }
    if (params.phix) {
        error += file_exists(params.phix, '--phix')
    }

    if (error > 0) {
        log.error('Cannot continue, please see --help for more information')
        exit 1
    }

    return fastq_type
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

def process_accessions(accession) {
    /* Parse line and determine if single end or paired reads*/
    return tuple(accession.trim(), "is_accession", [null, null])
}


def create_fastq_channel(fastq_input, fastq_type) {
    if (fastq_type == "fastqs") {
        return Channel.fromPath( file(fastq_input) )
            .splitCsv(header: true, sep: '\t')
            .map { row -> process_csv(row) }
    } else if (fastq_type == "ena_accessions") {
        return Channel.fromPath( file(params.accessions) )
            .splitText()
            .map { line -> process_accessions(line) }
    } else if (fastq_type == "ena_accession") {
        if (params.dry_run) {
            return [tuple("BACTOPIA_DRY_RUN", "is_accession", [null, null])]
        } else {
            return [tuple(params.accession, "is_accession", [null, null])]
        }
    } else if (fastq_type == "paired") {
        return [tuple(params.sample, false, [file(params.R1), file(params.R2)])]
    } else {
        return [tuple(params.sample, true, [file(params.SE)])]
    }
}

def create_reference_channel(channel_input) {
    return channel_input.map { it ->
        def split_references = []
        def sample = it[0]
        def single_end = it[1]
        def fastq = it[2]
        if (it[3] instanceof List) {
            it[3].each {
                split_references.add(tuple(sample, single_end, fastq, it))
            }
        } else {
            split_references.add(tuple(sample, single_end, fastq, it[3]))
        }
        return split_references
     }.flatMap { it -> L:{ it.collect { it } } }
}

def check_input_fastqs(fastq_input, fastq_type) {
    /* Read through --fastqs and verify each input exists. */
    if (fastq_type == "fastqs") {
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
                        log.error "LINE " + line + ':ERROR: Please verify ' + cols[1]+ ' exists, and try again'
                        error = true
                    }
                }
                if (cols[2]) {
                    if (!file(cols[2]).exists()) {
                        log.error "LINE " + line + ':ERROR: Please verify ' + cols[2]+ ' exists, and try again'
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
            log.error 'The header line (line 1) does not follow expected structure.'
        }

        if (error) {
            log.error 'Verify sample names are unique and/or FASTQ paths are correct'
            log.error 'See "--example_fastqs" for an example'
            log.error 'Exiting'
            exit 1
        }
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


def print_dataset_info(dataset_list, dataset_info) {
    log.info "Found ${dataset_list.size()} ${dataset_info}"
    dataset_list.each {
        log.info "\t${it}"
    }
}


def print_usage() {
    usage_text = params.help_all ? full_help() : basic_help()
    log.info"""
    ${PROGRAM_NAME} v${VERSION}
    ${basic_help()}
    ${params.help_all ? full_help() : ""}
    """.stripIndent()

    if (params.conda_help) {
        // Cleanup up the directory
        // This is only meant to be used with tests for conda build
        file("./work/").deleteDir()
        file("./.nextflow/").deleteDir()
        def files = new FileNameByRegexFinder().getFileNames('./', '.nextflow.log*')
        files.each { new File(it).delete()}
    }
    exit 0
}


def basic_help() {
    genome_size = params.genome_size ? params.genome_size : "Mash Estimate"
    return """
    Required Parameters:
        ### For Procesessing Multiple Samples
        --fastqs STR            An input file containing the sample name and
                                    absolute paths to FASTQs to process

        ### For Processing A Single Sample
        --R1 STR                First set of reads for paired end in compressed (gzip)
                                    FASTQ format

        --R2 STR                Second set of reads for paired end in compressed (gzip)
                                    FASTQ format

        --SE STR                Single end set of reads in compressed (gzip) FASTQ format

        --sample STR            The name of the input sequences

        ### For Downloading from ENA
        --accessions            An input file containing ENA/SRA experiement accessions to
                                    be processed

        --accession             A single ENA/SRA Experiment accession to be processed


    Dataset Parameters:
        --datasets DIR          The path to available datasets that have
                                    already been set up

        --species STR           Determines which species-specific dataset to
                                    use for the input sequencing

    Optional Parameters:
        --coverage INT          Reduce samples to a given coverage
                                    Default: ${params.coverage}x

        --genome_size INT       Expected genome size (bp) for all samples, a value of '0'
                                    will disable read error correction and read subsampling.
                                    Special values (requires --species):
                                        'min': uses minimum completed genome size of species
                                        'median': uses median completed genome size of species
                                        'mean': uses mean completed genome size of species
                                        'max': uses max completed genome size of species
                                    Default: Mash estimate

        --outdir DIR            Directory to write results to
                                    Default: ${params.outdir}

        --max_time INT          The maximum number of minutes a job should run before being halted.
                                    Default: ${params.max_time} minutes

        --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                    Default: ${params.max_memory} Gb

        --cpus INT              Number of processors made available to a single
                                    process.
                                    Default: ${params.cpus}

    Nextflow Related Parameters:
        --infodir DIR           Directory to write Nextflow summary files to
                                    Default: ${params.infodir}

        --condadir DIR          Directory to Nextflow should use for Conda environments
                                    Default: Bactopia's Nextflow directory

        --nfdir                 Print directory Nextflow has pulled Bactopia to

        --overwrite             Nextflow will overwrite existing output files.
                                    Default: ${params.overwrite}

        --conatainerPath        Path to Singularity containers to be used by the 'slurm'
                                    profile.
                                    Default: ${params.containerPath}

        --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                    will wait before execution.
                                    Default: ${params.sleep_time} seconds

    Useful Parameters:
        --available_datasets    Print a list of available datasets found based
                                    on location given by "--datasets"

        --example_fastqs        Print example of expected input for FASTQs file

        --check_fastqs          Verify "--fastqs" produces the expected inputs

        --compress              Compress (gzip) select outputs (e.g. annotation, variant calls)
                                    to reduce overall storage footprint.

        --keep_all_files        Keeps all analysis files created. By default, intermediate
                                    files are removed. This will not affect the ability
                                    to resume Nextflow runs, and only occurs at the end
                                    of the process.

        --dry_run               Mimics workflow execution, to help prevent errors realated to
                                    conda envs being built in parallel. Only useful on new
                                    installs of Bactopia.

        --version               Print workflow version information

        --help                  Show this message and exit

        --help_all              Show a complete list of adjustable parameters
    """
}


def full_help() {
    return """
    Additional Parameters:
    The description of the following parameters were taken from the program for
    which they apply to.

    Many of the default values were also taken from the program for which they
    apply to.

    ENA Download Parameters:
        --max_retry INT         Maximum times to retry downloads
                                    Default: ${params.max_retry}

        --use_ena               Download FASTQs from ENA with Aspera Connect.
                                    Default: Download from SRA

        --ftp_only              If "--use_ena" is enabled, FTP will be used to
                                    download FASTQs from ENA.

        --aspera_speed STR      Speed at which Aspera Connect will download.
                                    Default: ${params.aspera_speed}

    Estimate Genome Size Parameters:
        Only applied if the genome size is estimated.

        --min_genome_size INT   The minimum estimated genome size allowed for the input sequence
                                to continue downstream analyses.
                                Default: ${params.min_genome_size}

        --max_genome_size INT   The maximum estimated genome size allowed for the input sequence
                                to continue downstream analyses.
                                Default: ${params.max_genome_size}

    QC Reads Parameters:
        --qc_ram INT            Try to keep RAM usage below this many GB
                                    Default: ${params.qc_ram} GB

        --min_basepairs INT     The minimum amount of input sequenced basepairs required
                                    to continue downstream analyses.
                                    Default: ${params.min_basepairs}

        --adapters FASTA        Illumina adapters to remove
                                    Default: BBmap adapters

        --adapter_k INT         Kmer length used for finding adapters. Adapters
                                    shorter than k will not be found
                                    Default: ${params.adapter_k}

        --phix FASTA            phiX174 reference genome to remove
                                    Default: NC_001422

        --phix_k INT            Kmer length used for finding phiX174.
                                    Contaminants shorter than k will not be
                                    found
                                    Default: ${params.phix_k}

        --ktrim STR             Trim reads to remove bases matching reference
                                    kmers
                                    Values:
                                        f (do not trim)
                                        r (trim to the right, Default)
                                        l (trim to the left)

        --mink INT              Look for shorter kmers at read tips down to this
                                    length, when k-trimming or masking. 0 means
                                    disabled. Enabling this will disable
                                    maskmiddle
                                    Default: ${params.mink}

        --hdist INT             Maximum Hamming distance for ref kmers (subs only)
                                    Memory use is proportional to (3*K)^hdist
                                    Default: ${params.hdist}

        --tpe BOOL              When kmer right-trimming, trim both reads to the
                                    minimum length of either
                                    Values:
                                        f (do not equally trim)
                                        t (equally trim to the right, Default)

        --tbo BOOL              Trim adapters based on where paired reads overlap
                                    Values:
                                        f (do not trim by overlap)
                                        t (trim by overlap, Default)

        --qtrim STR             Trim read ends to remove bases with quality
                                    below trimq. Performed AFTER looking for
                                    kmers
                                    Values:
                                        rl (trim both ends, Default)
                                        f (neither end)
                                        r (right end only)
                                        l (left end only)
                                        w (sliding window)

        --trimq FLOAT           Regions with average quality BELOW this will be
                                    trimmed if qtrim is set to something other
                                    than f
                                    Default: ${params.trimq}

        --maq INT               Reads with average quality (after trimming)
                                    below this will be discarded
                                    Default: ${params.maq}

        --minlength INT         Reads shorter than this after trimming will be
                                    discarded. Pairs will be discarded if both
                                    are shorter
                                    Default: ${params.minlength}

        --ftm INT               If positive, right-trim length to be equal to
                                    zero, modulo this number
                                    Default: ${params.ftm}

        --tossjunk              Discard reads with invalid characters as bases
                                    Values:
                                        f (keep all reads)
                                        t (toss reads with ambiguous bases, Default)

        --qout STR              Output quality offset
                                    Values:
                                        33 (PHRED33 offset quality scores, Default)
                                        64 (PHRED64 offset quality scores)
                                        auto (keeps the current input offset)

        --xmx STR               This will be passed to Java to set memory usage
                                    Examples:
                                        '8g' will specify 8 gigs of RAM (Default)
                                        '20g' will specify 20 gigs of RAM
                                        '200m' will specify 200 megs of RAM

        --maxcor INT            Max number of corrections within a 20bp window
                                    Default: ${params.maxcor}

        --sampleseed INT        Set to a positive number to use as the rng seed
                                    for sampling
                                    Default: ${params.sampleseed}


    Assembly Parameters:
        --shovill_ram INT       Try to keep RAM usage below this many GB
                                    Default: ${params.shovill_ram} GB

        --assembler STR         Assembler: megahit velvet skesa spades
                                    Default: ${params.assembler}

        --min_contig_len INT    Minimum contig length <0=AUTO>
                                    Default: ${params.min_contig_len}

        --min_contig_cov INT    Minimum contig coverage <0=AUTO>
                                    Default: ${params.min_contig_cov}

        --contig_namefmt STR    Format of contig FASTA IDs in 'printf' style
                                    Default: ${params.contig_namefmt}

        --shovill_opts STR      Extra assembler options in quotes eg.
                                    spades: "--untrusted-contigs locus.fna" ...
                                    Default: ''

        --shovill_kmers STR     K-mers to use <blank=AUTO>
                                    Default: ''

        --trim                  Enable adaptor trimming

        --nostitch              Disable read stitching

        --nocorr                Disable post-assembly correction


    Count 31mers Parameters:
        --cortex_ram INT        Try to keep RAM usage below this many GB
                                    Default: ${params.cortex_ram} GB

        --keep_singletons       Keep all counted 31-mers
                                    Default: Filter out singletons

    Annotation Parameters:
        --centre STR            Sequencing centre ID
                                    Default: ''

        --addgenes              Add 'gene' features for each 'CDS' feature

        --addmrna               Add 'mRNA' features for each 'CDS' feature

        --rawproduct            Do not clean up /product annotation

        --cdsrnaolap            Allow [tr]RNA to overlap CDS

        --prokka_evalue STR     Similarity e-value cut-off
                                    Default: ${params.prokka_evalue}

        --prokka_coverage INT   Minimum coverage on query protein
                                     Default: ${params.prokka_coverage}

        --norrna                Don't run rRNA search

        --notrna                Don't run tRNA search

        --rnammer               Prefer RNAmmer over Barrnap for rRNA prediction

    Minmer Sketch Parameters:
        --mash_sketch INT       Sketch size. Each sketch will have at most this
                                    many non-redundant min-hashes.
                                    Default: ${params.mash_sketch}

        --sourmash_scale INT    Choose number of hashes as 1 in FRACTION of
                                    input k-mers
                                    Default: ${params.sourmash_scale}


    Minmer Query Parameters:
        --minmer_ram INT        Try to keep RAM usage below this many GB
                                    Default: ${params.minmer_ram} GB

        --screen_w              Winner-takes-all strategy for identity estimates.
                                    After counting hashes for each query, hashes
                                    that appear in multiple queries will be
                                    removed from all except the one with the best
                                    identity (ties broken by larger query), and
                                    other identities will be reduced. This
                                    removes output redundancy, providing a rough
                                    compositional outline.
                                    Default: True

        --screen_i FLOAT        Minimum identity to report. Inclusive unless set
                                    to zero, in which case only identities greater
                                    than zero (i.e. with at least one shared hash)
                                    will be reported. Set to -1 to output
                                    everything.
                                    Default: ${params.screen_i}

    Ariba Parameters:
        --nucmer_min_id INT     Minimum alignment identity (delta-filter -i)
                                    Default: ${params.nucmer_min_id}

        --nucmer_min_len INT    Minimum alignment length (delta-filter -i)
                                    Default: ${params.nucmer_min_len}

        --nucmer_breaklen INT   Value to use for -breaklen when running nucmer
                                    Default: ${params.nucmer_breaklen}

        --assembly_cov INT      Target read coverage when sampling reads for
                                    assembly
                                    Default: ${params.assembly_cov}

        --min_scaff_depth INT   Minimum number of read pairs needed as evidence
                                    for scaffold link between two contigs
                                    Default: ${params.min_scaff_depth}

        --spades_options STR    Extra options to pass to Spades assembler
                                    Default: ${params.spades_options}

        --assembled_threshold FLOAT (between 0 and 1)
                                If proportion of gene assembled (regardless of
                                    into how many contigs) is at least this
                                    value then the flag gene_assembled is set
                                    Default: ${params.assembled_threshold}

        --gene_nt_extend INT    Max number of nucleotides to extend ends of gene
                                    matches to look for start/stop codons
                                    Default: ${params.gene_nt_extend}

        --unique_threshold FLOAT (between 0 and 1)
                                If proportion of bases in gene assembled more
                                    than once is <= this value, then the flag
                                    unique_contig is set
                                    Default: ${params.unique_threshold}

        --ariba_no_clean        Do not clean up intermediate files created by
                                    Ariba. By default, the local assemblies are
                                    deleted.

    Call Variant Parameters:
        --snippy_ram INT        Try and keep RAM under this many GB
                                    Default: ${params.snippy_ram} GB

        --mapqual INT           Minimum read mapping quality to consider
                                    Default: ${params.mapqual}

        --basequal INT          Minimum base quality to consider
                                    Default: ${params.basequal}

        --mincov INT            Minimum site depth to for calling alleles
                                    Default: ${params.mincov}

        --minfrac FLOAT         Minimum proportion for variant evidence (0=AUTO)
                                    Default: ${params.minfrac}

        --minqual INT           Minimum QUALITY in VCF column 6
                                    Default: ${params.minqual}

        --maxsoft INT           Maximum soft clipping to allow
                                    Default: ${params.maxsoft}

        --bwaopt STR            Extra BWA MEM options, eg. -x pacbio
                                    Default: ''

        --fbopt STR             Extra Freebayes options,
                                    eg. --theta 1E-6 --read-snp-limit 2
                                    Default: ''

    Nearest Neighbor Reference Genomes:
        --max_references INT    Maximum number of nearest neighbor reference genomes to
                                    download for variant calling.
                                    Default: 1

        --random_tie_break      On references with matching distances, randomly select one.
                                    Default: Pick earliest accession number

        --disable_auto_variants Disable automatic selection of reference genome based on
                                    Mash distances.

    Insertion Sequence Parameters:
        --min_clip INT          Minimum size for softclipped region to be
                                    extracted from initial mapping
                                    Default: ${params.min_clip}

        --max_clip INT          Maximum size for softclipped regions to be
                                    included
                                    Default: ${params.max_clip}

        --cutoff INT            Minimum depth for mapped region to be kept in
                                    bed file
                                    Default: ${params.cutoff}

        --novel_gap_size INT    Distance in base pairs between left and right
                                    flanks to be called a novel hit
                                    Default: ${params.novel_gap_size}

        --min_range FLOAT       Minimum percent size of the gap to be called a
                                    known hit
                                    Default: ${params.min_range}

        --max_range FLOAT       Maximum percent size of the gap to be called a
                                    known hit
                                    Default: ${params.max_range}

        --merging INT           Value for merging left and right hits in bed
                                    files together to simply calculation of
                                    closest and intersecting regions
                                    Default: ${params.merging}

        --ismap_all             Switch on all alignment reporting for bwa

        --ismap_minqual INT     Mapping quality score for bwa
                                    Default: ${params.ismap_minqual}

    BLAST Parameters:
        --perc_identity INT     Percent identity
                                    Default: ${params.perc_identity}

        --qcov_hsp_perc INT     Percent query coverage per hsp
                                    Default: ${params.qcov_hsp_perc}

        --max_target_seqs INT   Maximum number of aligned sequences to
                                    keep
                                    Default: ${params.max_target_seqs}

        --outfmt STR            BLAST alignment view options
                                    Default: '${params.outfmt}'

    Mapping Parameters:
        --keep_unmapped_reads   Keep unmapped reads, this does not affect variant
                                    calling.

        --bwa_mem_opts STR      Extra BWA MEM options
                                    Default: ''

        --bwa_aln_opts STR      Extra BWA ALN options
                                    Default: ''

        --bwa_samse_opts STR    Extra BWA SAMSE options
                                    Default: ''

        --bwa_sampe_opts STR    Extra BWA SAMPE options
                                    Default: ''

        --bwa_n INT             Maximum number of alignments to output in the XA
                                    tag for reads paired properly. If a read has
                                    more than INT hits, the XA tag will not be
                                    written.
                                    Default: ${params.bwa_n}

    Antimicrobial Resistance Parameters:
        --update_amr            Force amrfinder to update its database.

        --amr_ident_min         Minimum identity for nucleotide hit (0..1). -1
                                    means use a curated threshold if it exists and
                                    0.9 otherwise
                                    Default: ${params.amr_ident_min}

        --amr_coverage_min      Minimum coverage of the reference protein (0..1)
                                    Default: ${params.amr_coverage_min}

        --amr_organism          Taxonomy group: Campylobacter, Escherichia, Klebsiella
                                    Salmonella, Staphylococcus, Vibrio
                                    Default: ''

        --amr_translation_table NCBI genetic code for translated BLAST
                                    Default: ${params.amr_translation_table}

        --amr_plus              Add the plus genes to the report

        --amr_report_common     Suppress proteins common to a taxonomy group
        
    """
}
