#! /usr/bin/env nextflow
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine
import java.nio.file.Path
import java.nio.file.Paths
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version
log.info "bactopia tool ${PROGRAM_NAME} - ${VERSION}"

// Validate parameters
if (params.help || workflow.commandLine.trim().endsWith(workflow.scriptName)) print_help();
if (params.version) print_version();
check_input_params()
samples = gather_sample_set(params.bactopia, params.exclude, params.sleep_time)

// Setup output directories
outdir = "${params.outdir}/bactopia-tool/phyloflash"

process reconstruct_16s {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "${sample}/*"
    
    input:
    set val(sample), val(single_end), file(fq) from Channel.fromList(samples)

    output:
    file "${sample}/*" 
    file "${sample}/${sample}.SSU.collection.fasta" into ALIGNMENT

    shell:
    if (single_end)
    """
    phyloFlash.pl -dbhome "!{params.phyloflash_db}" -read1 !{fq[0]} -read2 !{fq[1]} \
                  -lib !{sample} -CPUs !{task.cpus} -taxlevel !{params.taxlevel}
    mkdir !{sample}
    mv !{sample}.* !{sample}
    """
    else 
    """
    phyloFlash.pl -dbhome "!{params.phyloflash_db}" -read1 !{fq[0]} -read2 !{fq[1]} \
                  -lib !{sample} -CPUs !{task.cpus} -taxlevel !{params.taxlevel}
    mkdir !{sample}
    mv !{sample}.* !{sample}
    """
}

process align_16s {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "16s-alignment.fasta"
    
    input:
    file(fasta) from ALIGNMENT.collect()

    output:
    file "16s-alignment.fa" into TREE

    shell:
    """
    cat *.fasta > 16s-merged.fa
    mafft --thread !{task.cpus} !{params.mafft_opts} 16s-merged.fa > 16s-alignment.fa
    """
}

process create_phylogeny {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "iqtree/*"
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "${params.prefix}.iqtree"

    input:
    file fasta from TREE

    output:
    file 'iqtree/*'
    file "${params.prefix}.iqtree"

    shell:
    """
    mkdir iqtree
    iqtree -s !{fasta} -m !{params.m} -nt !{task.cpus} -pre iqtree/16s \
           -bb !{params.bb} -alrt !{params.alrt} -wbt -wbtl \
           -alninfo !{params.iqtree_opts}
    cp iqtree/16s.iqtree !{params.prefix}.iqtree
    """
}

workflow.onComplete {
    workDir = new File("${workflow.workDir}")
    workDirSize = toHumanString(workDir.directorySize())

    println """
    Bactopia Tool 'cgtree' - Execution Summary
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

    // Check for unexpected paramaters
    error += check_unknown_params()

    if (params.bactopia) {
        error += file_exists(params.bactopia, '--bactopia')
    } else if (params.exclude) {
        error += file_exists(params.exclude, '--exclude')
    } else {
        log.error """
        The required '--bactopia' and/or '--phyloflash_db' parameter is missing, please check and try again.

        Required Parameters:
            --bactopia STR          Directory containing Bactopia analysis results for all samples.

            --phyloflash_db STR     Directory containing a pre-built phyloFlash database.
        """.stripIndent()
        error += 1
    }

    error += is_positive_integer(params.cpus, 'cpus')
    error += is_positive_integer(params.max_time, 'max_time')
    error += is_positive_integer(params.max_memory, 'max_memory')
    error += is_positive_integer(params.sleep_time, 'sleep_time')
    error += is_positive_integer(params.taxlevel, 'taxlevel')
    error += is_positive_integer(params.bb, 'bb')
    error += is_positive_integer(params.alrt, 'alrt')

    if (error > 0) {
        log.error('Cannot continue, please see --help for more information')
        exit 1
    }
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

def is_sample_dir(sample, dir){
    return file("${dir}/${sample}/${sample}-genome-size.txt").exists()
}

def build_fastq_tuple(sample, dir) {
    se = "${dir}/${sample}/quality-control/${sample}.fastq.gz"
    pe1 = "${dir}/${sample}/quality-control/${sample}_R1.fastq.gz"
    pe2 = "${dir}/${sample}/quality-control/${sample}_R2.fastq.gz"
    if (file(se).exists()) {
        return tuple(sample, true, [file(se)])
    } else if (file(pe1).exists() && file(pe2).exists()) {
        return tuple(sample, false, [file(pe1), file(pe2)])
    } else {
        log.error("Could not locate FASTQs for ${sample}, please verify existence. Unable to continue.")
        exit 1
    }    
}

def gather_sample_set(bactopia_dir, exclude_list, sleep_time) {
    exclusions = []
    if (exclude_list) {
        new File(exclude_list).eachLine { line ->
            exclusions << line.trim()
        }
        log.info "Excluding ${exclusions.size} samples from the analysis"
    }


    sample_list = []
    file(bactopia_dir).eachFile { item ->
        
        if( item.isDirectory() ) {
            sample = item.getName()
            if (!exclusions.contains(sample)) {
                if (is_sample_dir(sample, bactopia_dir)) {
                    sample_list << build_fastq_tuple(sample, bactopia_dir)
                }
            }
        }
    }
    log.info "Found ${sample_list.size} samples to process"
    log.info "\nIf this looks wrong, now's your chance to back out (CTRL+C 3 times)."
    log.info "Sleeping for ${sleep_time} seconds..."
    sleep(sleep_time * 1000)
    return sample_list[1..5]
}


def print_help() {
    log.info"""
    Required Parameters:
        --bactopia STR          Directory containing Bactopia analysis results for all samples.

        --phyloflash_db STR     Directory containing a pre-built phyloFlash database.

    Optional Parameters:
        --exclude STR           A text file containing sample names to exclude from the
                                    pangenome. The expected format is a single sample per line.

        --prefix DIR            Prefix to use for final output files
                                    Default: ${params.prefix}

        --outdir DIR            Directory to write results to
                                    Default: ${params.outdir}

        --max_time INT          The maximum number of minutes a job should run before being halted.
                                    Default: ${params.max_time} minutes

        --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                    Default: ${params.max_memory} Gb

        --cpus INT              Number of processors made available to a single
                                    process.
                                    Default: ${params.cpus}

    phyloFlash Related Parameters:
        --taxlevel INT          Level in the taxonomy string to summarize read counts per taxon.
                                    Numeric and 1-based (i.e. "1" corresponds to "Domain").
                                    Default: ${params.taxlevel}

        --phyloflash_opts STR   Extra phyloFlash options in quotes.
                                    Default: ''


    MAFFT Related Parameters:
        --mafft_opts STR        MAFFT options to include (in quotes).
                                    Default: ''

    IQ-TREE Related Parameters:
        --skip_phylogeny        Skip the creation a core-genome based phylogeny

        --m STR                 Substitution model name
                                    Default: ${params.m}

        --bb INT                Ultrafast bootstrap replicates
                                    Default: ${params.bb}

        --alrt INT              SH-like approximate likelihood ratio test replicates
                                    Default: ${params.alrt}

        --asr                   Ancestral state reconstruction by empirical Bayes
                                    Default: ${params.asr}

        --iqtree_opts STR       Extra IQ-TREE options in quotes.
                                    Default: ''

    Nextflow Related Parameters:
        --infodir DIR           Directory to write Nextflow summary files to
                                    Default: ${params.infodir}

        --overwrite             Nextflow will overwrite existing output files.
                                    Default: ${params.overwrite}

        --conatainerPath        Path to Singularity containers to be used by the 'slurm'
                                    profile.
                                    Default: ${params.containerPath}

        --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                    will wait before execution.
                                    Default: ${params.sleep_time} seconds
    Useful Parameters:
        --version               Print workflow version information
        --help                  Show this message and exit
    """.stripIndent()
    exit 0
}
