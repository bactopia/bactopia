#! /usr/bin/env nextflow
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version
OUTDIR = "${params.outdir}/bactopia-tools/${PROGRAM_NAME}/${params.prefix}"
OVERWRITE = workflow.resume || params.force ? true : false

// Validate parameters
if (params.version) print_version();
log.info "bactopia tools ${PROGRAM_NAME} - ${VERSION}"
if (params.help || workflow.commandLine.trim().endsWith(workflow.scriptName)) print_help();
check_input_params()
samples = gather_sample_set(params.bactopia, params.exclude, params.include, params.sleep_time, params.insertions)

process collect_reference {    
    publishDir "${OUTDIR}", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "genbank/${reference_name}.gbk"

    input:
    file(reference) from Channel.fromPath(params.reference).ifEmpty { "EMPTY.gbk" }

    output:
    file "genbank/${reference_name}.gbk" optional true
    file("genbank/${reference_name}.gbk") into ISMAPPER_REFERENCE

    shell:
    is_compressed = null
    section = null
    reference_name = null
    if (reference) {
        is_compressed = params.reference.endsWith(".gz") ? true : false
        reference_name = reference.getSimpleName()
    } else {
        section = params.accession.startsWith('GCF') ? 'refseq' : 'genbank'
        reference_name = params.accession.split(/\./)[0]
    }
    """
    mkdir genbank
    if [ "!{params.accession}" != "null" ]; then
        ncbi-genome-download bacteria -l complete -o ./ -F genbank \
                                      -s !{section} -A !{params.accession} -r 50

        find . -name "*!{params.accession}*.gbff.gz" | xargs -I {} mv {} genbank/
        if [ "!{section}" == 'refseq' ]; then
            rename 's/(GCF_\\d+).*/\$1.gbk.gz/' genbank/*
        else
            rename 's/(GCA_\\d+).*/\$1.gbk.gz/' genbank/*
        fi
        gunzip genbank/*
    else
        if [ "!{is_compressed}" == "true" ]; then
            zcat !{reference} > genbank/!{reference_name}.gbk
        else 
            cat !{reference} > genbank/!{reference_name}.gbk
        fi
    fi
    """
}

process insertion_sites {
    /* Query a set of insertion sequences (FASTA) against annotated GenBank file using ISMapper.  */
    tag "${sample} - ${insertion_name}"
    publishDir "${OUTDIR}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${sample}/*"

    input:
    set val(sample), file(fq), file(insertions) from Channel.fromList(samples)
    each file(reference) from ISMAPPER_REFERENCE

    output:
    file("${sample}/*")

    shell:
    insertion_name = insertions.getSimpleName()
    all = params.ismap_all ? "--a" : ""
    """
    ismap --reads !{sample}_R*.fastq.gz \
        --queries !{insertions} \
        --reference !{reference} \
        --log !{sample}-!{insertion_name} \
        --min_clip !{params.min_clip} \
        --max_clip !{params.max_clip} \
        --cutoff !{params.cutoff} \
        --novel_gap_size !{params.novel_gap_size} \
        --min_range !{params.min_range} \
        --max_range !{params.max_range} \
        --merging !{params.merging} \
        --T !{params.ismap_minqual} \
        --t !{task.cpus} !{all}

    mv !{sample}-!{insertion_name}.log !{sample}/
    """
}


workflow.onComplete {
    workDir = new File("${workflow.workDir}")
    workDirSize = toHumanString(workDir.directorySize())

    println """
    Bactopia Tool '${PROGRAM_NAME}' - Execution Summary
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
    log.info "bactopia tools ${PROGRAM_NAME} - ${VERSION}"
    exit 0
}

def file_exists(file_name, parameter) {
    if (!file(file_name).exists()) {
        log.error('Invalid input ('+ parameter +'), please verify "' + file_name + '" exists.')
        return 1
    }
    return 0
}

def output_exists(outdir, force, resume) {
    if (!resume && !force) {
        if (file(OUTDIR).exists()) {
            files = file(OUTDIR).list()
            total_files = files.size()
            if (total_files == 1) {
                if (files[0] != 'bactopia-info') {
                    return 1
                }
            } else if (total_files > 1){
                return 1
            }
        }
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

    IGNORE_LIST = ['container-path']
    params.each { k,v ->
        if (!valid_params.contains(k)) {
            if (!IGNORE_LIST.contains(k)) {
                log.error("'--${k}' is not a known parameter")
                error = 1
            }
        }
    }

    return error
}

def check_input_params() {
    // Check for unexpected paramaters
    error = check_unknown_params()

    if (params.bactopia && params.insertions && params.reference) {
        error += file_exists(params.bactopia, '--bactopia')
        error += file_exists(params.insertions, '--insertions')
        error += file_exists(params.reference, '--reference')
    } else if (params.bactopia && params.insertions && params.accession) {
        error += file_exists(params.bactopia, '--bactopia')
        error += file_exists(params.insertions, '--insertions')
    } else {
        log.error """
        Missing one or more required parameters, please check and try again.

        Required Parameters:
            --bactopia STR          Directory containing Bactopia analysis results for all samples.

            --insertions STR        Multifasta file with insertion sequence(s) to be mapped to.

            *Note: --reference or --accession is required.*
            --reference STR         Reference genome for typing against in GenBank format.

            --accession STR         The Assembly accession (e.g. GC(A|F)*.*) of the reference to
                                        download from RefSeq.
        """.stripIndent()
        error += 1
    }

    if (params.include) {
        error += file_exists(params.include, '--include')
    } 

    if (params.exclude) {
        error += file_exists(params.exclude, '--exclude')
    } 

    error += is_positive_integer(params.cpus, 'cpus')
    error += is_positive_integer(params.max_time, 'max_time')
    error += is_positive_integer(params.max_memory, 'max_memory')
    error += is_positive_integer(params.sleep_time, 'sleep_time')

    error += is_positive_integer(params.min_clip, 'min_clip')
    error += is_positive_integer(params.max_clip, 'max_clip')
    error += is_positive_integer(params.cutoff, 'cutoff')
    error += is_positive_integer(params.novel_gap_size, 'novel_gap_size')
    error += is_positive_integer(params.min_range, 'min_range')
    error += is_positive_integer(params.max_range, 'max_range')
    error += is_positive_integer(params.merging, 'merging')
    error += is_positive_integer(params.ismap_minqual, 'ismap_minqual')

    // Check for existing output directory
    if (output_exists(OUTDIR, params.force, workflow.resume)) {
        log.error("Output directory (${OUTDIR}) exists, Bactopia will not continue unless '--force' is used.")
        error += 1
    }

    // Check publish_mode
    ALLOWED_MODES = ['copy', 'copyNoFollow', 'link', 'rellink', 'symlink']
    if (!ALLOWED_MODES.contains(params.publish_mode)) {
        log.error("'${params.publish_mode}' is not a valid publish mode. Allowed modes are: ${ALLOWED_MODES}")
        error += 1
    }

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
        if (!value.toString().isNumber()) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not numeric.')
            error = 1
        } else if (value.toString().toFloat() < 0) {
            log.error('Invalid input (--'+ name +'), "' + value + '"" is not positive.')
            error = 1
        }
    }
    return error
}

def is_sample_dir(sample, dir){
    return file("${dir}/${sample}/${sample}-genome-size.txt").exists()
}

def build_fastq_tuple(sample, dir, insertions) {
    se = "${dir}/${sample}/quality-control/${sample}.fastq.gz"
    pe1 = "${dir}/${sample}/quality-control/${sample}_R1.fastq.gz"
    pe2 = "${dir}/${sample}/quality-control/${sample}_R2.fastq.gz"

    if (file(se).exists()) {
        log.info "Excluding ${sample} from the analysis, must be paired-end"
    } else if (file(pe1).exists() && file(pe2).exists()) {
        return tuple(sample, [file(pe1), file(pe2)], file(insertions))
    } else {
        log.error("Could not locate FASTQs for ${sample}, please verify existence. Unable to continue.")
        exit 1
    }
}

def gather_sample_set(bactopia_dir, exclude_list, include_list, sleep_time, insertions) {
    include_all = true
    inclusions = []
    exclusions = []
    IGNORE_LIST = ['.nextflow', 'bactopia-info', 'bactopia-tools', 'work',]
    if (include_list) {
        new File(include_list).eachLine { line -> 
            inclusions << line.trim()
        }
        include_all = false
        log.info "Including ${inclusions.size} samples for analysis"
    }
    else if (exclude_list) {
        new File(exclude_list).eachLine { line -> 
            exclusions << line.trim().split('\t')[0]
        }
        log.info "Excluding ${exclusions.size} samples from the analysis"
    }
    
    sample_list = []
    file(bactopia_dir).eachFile { item ->
        if( item.isDirectory() ) {
            sample = item.getName()
            if (!IGNORE_LIST.contains(sample)) {
                if (inclusions.contains(sample) || include_all) {
                    if (!exclusions.contains(sample)) {
                        if (is_sample_dir(sample, bactopia_dir)) {
                            sample_list << build_fastq_tuple(sample, bactopia_dir, insertions)
                        } else {
                            log.info "${sample} is missing genome size estimate file"
                        }
                    }
                }
            }
        }
    }

    log.info "Found ${sample_list.size} samples to process"
    log.info "\nIf this looks wrong, now's your chance to back out (CTRL+C 3 times)."
    log.info "Sleeping for ${sleep_time} seconds..."
    sleep(sleep_time * 1000)
    return sample_list
}


def print_help() {
    log.info"""
    Required Parameters:
        --bactopia STR          Directory containing Bactopia analysis results for all samples.

        --insertions STR        Multifasta file with insertion sequence(s) to be mapped to.

        *Note: --reference or --accession is required.*
        --reference STR         Reference genome for typing against in GenBank format.

        --accession STR         The Assembly accession (e.g. GC(A|F)*.*) of the reference to
                                    download from RefSeq.

    Optional Parameters:
        --include STR           A text file containing sample names to include in the
                                    analysis. The expected format is a single sample per line.

        --exclude STR           A text file containing sample names to exclude from the
                                    analysis. The expected format is a single sample per line.

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

    ISMapper Parameters:
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

    Nextflow Related Parameters:
        --condadir DIR          Directory to Nextflow should use for Conda environments
                                    Default: Bactopia's Nextflow directory
        
        --publish_mode          Set Nextflow's method for publishing output files. Allowed methods are:
                                    'copy' (default)    Copies the output files into the published directory.

                                    'copyNoFollow' Copies the output files into the published directory 
                                                   without following symlinks ie. copies the links themselves.

                                    'link'    Creates a hard link in the published directory for each 
                                              process output file.

                                    'rellink' Creates a relative symbolic link in the published directory
                                              for each process output file.

                                    'symlink' Creates an absolute symbolic link in the published directory 
                                              for each process output file.

                                    Default: ${params.publish_mode}

        --force                 Nextflow will overwrite existing output files.
                                    Default: ${params.force}

        --conatainerPath        Path to Singularity containers to be used by the 'slurm'
                                    profile.
                                    Default: ${params.containerPath}

        --sleep_time            After reading datases, the amount of time (seconds) Nextflow
                                    will wait before execution.
                                    Default: ${params.sleep_time} seconds

        --nfconfig STR          A Nextflow compatible config file for custom profiles. This allows 
                                    you to create profiles specific to your environment (e.g. SGE,
                                    AWS, SLURM, etc...). This config file is loaded last and will 
                                    overwrite existing variables if set.
                                    Default: Bactopia's default configs

        -resume                 Nextflow will attempt to resume a previous run. Please notice it is 
                                    only a single '-'

    Useful Parameters:
        --version               Print workflow version information
        --help                  Show this message and exit
    """.stripIndent()
    exit 0
}
