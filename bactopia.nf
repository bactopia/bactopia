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

// Set the maximum number of cpus to use
config.poolSize = params.max_cpus.toInteger()
cpus = params.cpus
if (cpus > params.max_cpus) {
    cpus = params.max_cpus
    log.info "--cpus ${params.cpus} exceeded --max_cpus ${params.max_cpus}, changed ${params.cpus} to ${cpus}"
}

// Setup output directories
outdir = params.outdir ? params.outdir : './'

process qc_reads {
    /* Cleanup the reads using Illumina-Cleanup */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from create_fastq_channel(params.fastqs)

    output:
    file "quality-control/*"
    set val(sample), val(single_end), file("quality-control/${sample}*.fastq.gz") into ASSEMBLY

    shell:
    template(task.ext.template)
}

process assemble_genome {
    /* Assemble the genome using Shovill, SKESA is used by default */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}/assembly", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from ASSEMBLY

    output:
    file "shovill*"
    file "${sample}.fna.gz"
    file "${sample}.fna.json"
    set val(sample), file("${sample}.fna.gz") into ANNOTATION

    shell:
    template(task.ext.template)
}

process annotate_genome {
    /* Annotate the assembly using Prokka, use a proteins fasta if available */
    cpus cpus
    tag "${sample}"
    publishDir "${outdir}/${sample}", mode: 'copy', overwrite: true

    input:
    set val(sample), file(fasta) from ANNOTATION

    output:
    file 'annotation/*'

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    proteins = prokka_proteins ? "--proteins !{prokka_proteins}" : ""
    template(task.ext.template)
}


process count_31mers {
    /* Count 31mers in the reads using McCortex */
    cpus cpus
    tag "${sample} - ${database}"
    publishDir "${outdir}/${sample}/kmers", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from COUNT_31MERS

    output:
    file "${sample}.ctx"

    shell:
    template(task.ext.template)

}

process ariba_databases {
    /* Run reads against all available (if any) Ariba databases */
    cpus cpus
    tag "${sample} - ${database}"
    publishDir "${outdir}/${sample}/ariba", mode: 'copy', overwrite: true

    input:
    set val(sample), val(single_end), file(fq) from ARIBA_DATABASES
    each database from available_ariba_databases

    output:
    file "${database}/*"

    shell:
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

def clean_cache() {
    // No need to resume completed run so remove cache.
    file('./work/').deleteDir()
    file('./.nextflow/').deleteDir()
    file('.nextflow.log').delete()
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

    if (!is_positive_integer(params.genome_size, 'genome_size')) {
        error = true
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
