#! /usr/bin/env nextflow
import java.nio.file.Paths
BACTOPIA_DIR = params.bactopia
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version
OUTDIR = "${params.outdir}/bactopia-tools/${PROGRAM_NAME}/${params.prefix}"
OVERWRITE = workflow.resume || params.force ? true : false

// Validate parameters
if (params.version) print_version();
log.info "bactopia tools ${PROGRAM_NAME} - ${VERSION}"
if (params.help || workflow.commandLine.trim().endsWith(workflow.scriptName)) print_help();
check_input_params()

process summary {
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "${params.prefix}*"

    output:
    file "${params.prefix}*"
    file "${params.prefix}-exclude.txt" into ARIBA, AMRFINDER

    when:

    shell:
    verbose = params.verbose? "--verbose" : ""
    min_genome_size = params.min_genome_size? "--min_genome_size ${params.min_genome_size}" : ""
    max_genome_size = params.max_genome_size? "--max_genome_size ${params.max_genome_size}" : ""
    """
    bactopia-summary.py !{BACTOPIA_DIR} !{min_genome_size} !{max_genome_size} !{verbose} \
        --prefix !{params.prefix} \
        --gold_coverage !{params.gold_coverage} \
        --gold_quality !{params.gold_quality} \
        --gold_read_length !{params.gold_read_length} \
        --gold_contigs !{params.gold_contigs} \
        --silver_coverage !{params.silver_coverage} \
        --silver_quality !{params.silver_quality} \
        --silver_read_length !{params.silver_read_length} \
        --silver_contigs !{params.silver_contigs} \
        --min_coverage !{params.min_coverage} \
        --min_quality !{params.min_quality} \
        --min_read_length !{params.min_read_length} \
        --max_contigs !{params.max_contigs}
    """
}

process ariba_summary {
    publishDir "${OUTDIR}/ariba", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "ariba-*-summary.txt"

    input:
    file(exclude) from ARIBA

    output:
    file "ariba-*-summary.txt" optional true

    shell:
    all_hits = params.all_hits ? "--include_all" : ""
    verbose = params.verbose? "--verbose" : ""
    """
    ariba-summary.py !{BACTOPIA_DIR} --exclude !{exclude} !{all_hits} !{verbose}
    """
}

process amrfinder_summary {
    publishDir "${OUTDIR}/amrfinder", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "amrfinder-*-summary.txt"

    input:
    file(exclude) from AMRFINDER

    output:
    file "amrfinder-*-summary.txt" optional true

    shell:
    subclass = params.subclass ? "--subclass" : ""
    verbose = params.verbose? "--verbose" : ""
    """
    amrfinder-summary.py !{BACTOPIA_DIR} --exclude !{exclude} !{subclass} !{verbose}
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
    // Check for unexpected paramaters
    error = check_unknown_params()

    if (params.bactopia) {
        error += file_exists(params.bactopia, '--bactopia')
    } else {
        log.error """
        The required '--bactopia' parameter is missing, please check and try again.

        Required Parameters:
            --bactopia STR          Directory containing Bactopia analysis results for all samples.
        """.stripIndent()
        error += 1
    }

    error += is_positive_integer(params.cpus, 'cpus')
    error += is_positive_integer(params.max_time, 'max_time')
    error += is_positive_integer(params.max_memory, 'max_memory')
    error += is_positive_integer(params.sleep_time, 'sleep_time')

    error += is_positive_integer(params.gold_coverage, 'gold_coverage')
    error += is_positive_integer(params.gold_quality, 'gold_quality')
    error += is_positive_integer(params.gold_read_length, 'gold_read_length')
    error += is_positive_integer(params.gold_contigs, 'gold_contigs')
    error += is_positive_integer(params.silver_coverage, 'silver_coverage')
    error += is_positive_integer(params.silver_quality, 'silver_quality')
    error += is_positive_integer(params.silver_read_length, 'silver_read_length')
    error += is_positive_integer(params.silver_contigs, 'silver_contigs')
    error += is_positive_integer(params.min_coverage, 'min_coverage')
    error += is_positive_integer(params.min_quality, 'min_quality')
    error += is_positive_integer(params.min_read_length, 'min_read_length')
    error += is_positive_integer(params.max_contigs, 'max_contigs')

    if (params.min_genome_size) {
        error += is_positive_integer(params.min_genome_size, 'min_genome_size')
    }

    if (params.max_genome_size) {
        error += is_positive_integer(params.max_genome_size, 'max_genome_size')
    }
    
    Path bactopia_path = Paths.get(BACTOPIA_DIR); 
    if (!bactopia_path.isAbsolute()) {
        BACTOPIA_DIR = "${workflow.launchDir}/${BACTOPIA_DIR}"
        log.info("A relative path to the bactopia directory (--bactopia ${params.bactopia}) was given. Using ${BACTOPIA_DIR} instead.")
    }

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

def print_help() {
    log.info"""
    Required Parameters:
        --bactopia STR          Directory containing Bactopia analysis results for all samples.

    Bactopia Summary Parameters:
        --gold_coverage FLOAT   Minimum amount of coverage required for Gold status
                                    Default: ${params.gold_coverage}

        --gold_quality INT      Minimum per-read mean quality score required for Gold
                                    status 
                                    Default: ${params.gold_quality}

        --gold_read_length INT  Minimum mean read length required for Gold status
                                    Default: ${params.gold_read_length}

        --gold_contigs INT      Maximum contig count required for Gold status
                                    Default: ${params.gold_contigs}

        --silver_coverage FLOAT Minimum amount of coverage required for Silver status
                                    Default: ${params.silver_coverage}

        --silver_quality INT    Minimum per-read mean quality score required for
                                    Silver status 
                                    Default: ${params.silver_quality}

        --silver_read_length INT
                                Minimum mean read length required for Silver status
                                    Default: ${params.silver_read_length}

        --silver_contigs INT    Maximum contig count required for Silver status
                                    Default: ${params.silver_contigs}

        --min_coverage FLOAT    Minimum amount of coverage required to pass 
                                    Default: ${params.min_coverage}

        --min_quality INT       Minimum per-read mean quality score required to pass
                                    Default: ${params.min_quality}

        --min_read_length INT   Minimum mean read length required to pass 
                                    Default: ${params.min_read_length}

        --max_contigs INT       Maximum contig count required to pass
                                    Default: ${params.max_contigs}

        --min_genome_size INT   Minimum assembled genome size.
                                    Default: ${params.min_genome_size}

        --max_genome_size INT   Maximum assembled genome size.
                                    Default: ${params.max_genome_size}

    Ariba Summary Parameters:
        --all_hits              Include all hits (matches and partials) in the summary
                                    Default: Only report hits that are a match

    AMRFinder+ Summary Parameters:
        --subclass              Group the report by subclass (ex. Streptomycin).
                                    Default: Group by class (ex. Aminoglycoside)

    Optional Parameters:
        --prefix STR            Prefix to use for final output files
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

    Nextflow Related Parameters:
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
        --verbose               Increase the verbosity of processes.
        --version               Print workflow version information
        --help                  Show this message and exit
    """.stripIndent()
    exit 0
}
