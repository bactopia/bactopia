#! /usr/bin/env nextflow
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine
import java.nio.file.Path
import java.nio.file.Paths
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version
log.info "bactopia tool ${PROGRAM_NAME} - ${VERSION}"

// Validate parameters
if (params.help) print_usage();
if (workflow.commandLine.trim().endsWith(workflow.scriptName)) print_usage();
if (params.version) print_version();
check_input_params()

// Setup output directories
outdir = params.outdir ? params.outdir : './'
Channel
    .fromPath("${params.bactopia}/**/annotation/*.gff*")
    .count()
    .subscribe { println "Found ${it} GFF files for analyis" }
sleep(2000)
log.info "\nIf this looks wrong, now's your chance to back out (CTRL+C 3 times)."
log.info "Sleeping for ${params.sleep_time} seconds..."
sleep(params.sleep_time * 1000)

process pangenome {
    cpus !{params.cpus}
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "roary/*"

    input:
    file(gff) from Channel.fromPath("${params.bactopia}/**/annotation/*.gff*").collect()

    output:
    file 'roary/*'
    file 'roary/core_gene_alignment.aln' into START_TREE

    shell:
    n = params.n ? "-n" : ""
    s = params.s ? "-s" : ""
    ap = params.ap ? "-ap" : ""
    compressed = params.compressed ? "gunzip -f *.gff.gz" : "echo uncompressed"
    """
    !{compressed}
    roary -f roary -e !{n} -v -p !{task.cpus} !{s} !{ap} -g !{params.g} \
          -i !{params.i} -cd !{params.cd} -iv !{params.iv} -r *.gff
    """
}

process start_tree {
    cpus !{params.cpus}
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "start-tree/*"

    input:
    file fasta from START_TREE

    output:
    file 'start-tree/*'
    set file('start-tree/start-tree.treefile'), file(fasta) into RECOMBINATION

    shell:
    """
    mkdir start-tree
    iqtree -s !{fasta} -m !{params.m} -nt !{task.cpus} -fast -pre start-tree/start-tree
    """
}

process recombination {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "clonalframe/*"

    input:
    set file(start_tree), file(fasta) from RECOMBINATION

    output:
    file 'clonalframe/*'
    file 'clonalframe/core_gene_alignment-cfmasked.aln' into FINAL_TREE, SNP_DISTS

    shell:
    if (params.skip_clonalframe)
    """
    mkdir clonalframe
    cp !{fasta} clonalframe/core_gene_alignment-cfmasked.aln
    touch clonalframe/clonalframe-was-skipped.txt
    """
    else
    """
    mkdir clonalframe
    ClonalFrameML !{start_tree} !{fasta} clonalframe/clonalframe \
        -emsim !{params.emsim} !{params.clonal_opts}

    maskrc-svg.py clonalframe/clonalframe --aln !{fasta} --symbol '-' \
        --out clonalframe/core_gene_alignment-cfmasked.aln
    """
}

process final_tree {
    cpus !{params.cpus}
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "final-tree/*"

    input:
    file fasta from FINAL_TREE

    output:
    file 'final-tree/*'

    shell:
    """
    mkdir final-tree
    iqtree -s !{fasta} -m !{params.m} -nt !{task.cpus} -pre final-tree/final-tree \
           -bb !{params.bb} -alrt !{params.alrt} -wbt -wbtl \
           -alninfo !{params.iqtree_opts}
    """
}

process snp_dists {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite

    input:
    file fasta from SNP_DISTS

    output:
    file 'pairwise-snp-distance.txt'

    shell:
    b = params.b ? "" : "-b"
    """
    snp-dists !{b} !{fasta} > pairwise-snp-distance.txt
    """
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
    error += is_positive_integer(params.i, 'snippy_ram')
    error += is_positive_integer(params.cd, 'cortex_ram')
    error += is_positive_integer(params.g, 'qc_ram')
    error += is_positive_integer(params.emsim, 'minmer_ram')
    error += is_positive_integer(params.bb, 'snippy_ram')
    error += is_positive_integer(params.alrt, 'cortex_ram')

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

def help() {
    return """
    Required Parameters:
        --bactopia STR          Directory containing Bactopia analysis results for all samples.

    Optional Parameters:
        --outdir DIR            Directory to write results to
                                    Default: ${params.outdir}

        --max_time INT          The maximum number of minutes a job should run before being halted.
                                    Default: ${params.max_time} minutes

        --max_memory INT        The maximum amount of memory (Gb) allowed to a single process.
                                    Default: ${params.max_memory} Gb

        --cpus INT              Number of processors made available to a single
                                    process.
                                    Default: ${params.cpus}

        --compressed            Input GFFs are compressed (gzip)
                                    Default: ${params.compressed}

    Roary Related Parameters:
        --o STR                 Clusters output filename
                                    Default: ${params.o}

        --n                     Execute a fast core gene alignment with MAFFT
                                    Default: Use PRANK

        --i INT                 Minimum percentage identity for blastp
                                    Default: ${params.i}

        --cd INT                Percentage of isolates a gene must be in to be core
                                    Default: ${params.cd}%

        --g INT                 Maximum number of clusters
                                    Default: ${params.g}

        --s                     Do not split paralogs
                                    Default: ${params.s}

        --ap                    Allow paralogs in core alignment
                                    Default: ${params.ap}

        --iv STR                Change the MCL inflation value
                                    Default: ${params.iv}

    IQ-TREE Related Parameters:
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

    ClonalFrameML Related Parameters:
        --skip_clonalframe      Skip the ClonalFrameML and use the original core-genome
                                    alignment for the final tree.
                                    Default: ${params.skip_clonalframe}

        --emsim INT             Number of simulations to estimate uncertainty in the EM results.
                                    Default: ${params.emsim}

        --clonal_opts STR       Extra ClonalFrameML options in quotes.
                                    Default: ''

    SNP-Dists Related Parameters:
        --a                     Count all differences not just [AGTC]
                                    Default: ${params.a}

        --b                     Blank top left corner cell
                                    Default: ${params.b}

        --c                     Output CSV instead of TSV
                                    Default: ${params.c}

        --k                     Keep case, don't uppercase all letters
                                    Default: ${params.k}

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
    """
}
