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
is_compressed = check_gffs_exist(params.bactopia, params.sleep_time)

// Setup output directories
outdir = "${params.outdir}/bactopia-tool/roary"

process download_references {
    publishDir "${outdir}/refseq", mode: 'copy', overwrite: params.overwrite, pattern: "fasta/*.fna"
    
    input:
    file(gff) from create_gff_channel(params.bactopia, true).collect()
    val(is_compressed) from is_compressed

    output:
    file 'fasta/*.fna' into ANNOTATE

    when:
    species != null

    shell:
    """
    ncbi-genome-download bacteria -l complete -o ./ -F fasta -p !{task.cpus} \
                                  --genus "!{params.species}" -r 50
    mkdir fasta
    find -name "*.fna.gz" | xargs -I {} mv {} fasta/
    rename 's/(GCF_\d+).*/$1.fna.gz/' fasta/*
    gunzip fasta/*
    """

}

process annotate_references {
    publishDir "${outdir}/refseq/gff", mode: 'copy', overwrite: params.overwrite, pattern: "*.gff"

    input:
    file fasta from ANNOTATE

    output:
    file '*.gff' into REFERENCE_GFF

    shell:
    """
    prokka !{fasta} --cpus !{task.cpus} \
        --evalue '!{params.prokka_evalue}' \
        --coverage !{params.prokka_coverage}
    """

}

process build_pangenome {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "roary/*"
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "excluded_gff/*"

    input:
    file(sample_gff) from create_gff_channel(params.bactopia, true).collect()
    file(reference_gff) from REFERENCE_GFF.collect()
    file(exclude) from params.exclude
    val(is_compressed) from is_compressed

    output:
    file 'roary/*'
    file 'excluded_gff/*'
    file 'alignment.fa' into RECOMBINATION

    shell:
    n = params.n ? "-n" : ""
    s = params.s ? "-s" : ""
    ap = params.ap ? "-ap" : ""
    gunzip = is_compressed ? "gunzip -f *.gff.gz" : "echo uncompressed"
    filter = exclude.name != 'NO_FILE' ? "${exclude.name}" : '/dev/null'
    """
    !{gunzip}
    mkdir excluded_gff
    ls *.gff | grep -f !{filter} | xargs -I {} mv {} excluded_gff
    roary -f roary -e !{n} -v -p !{task.cpus} !{s} !{ap} -g !{params.g} \
          -i !{params.i} -cd !{params.cd} -iv !{params.iv} -r *.gff
    cp roary/core_gene_alignment.aln alignment.fa
    pigz -n --best -p !{task.cpus} roary/core_gene_alignment.aln
    """
}

process identify_recombination {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "clonalframe/*"
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "${params.prefix}.aligned.fa.gz"

    input:
    file fasta from RECOMBINATION

    output:
    file 'clonalframe/*' optional true
    file "${params.prefix}.aligned.fa.gz"
    file 'alignment-masked.fa' into FINAL_TREE, SNP_DISTS

    when:
    params.skip_phylogeny == false

    shell:
    if (params.skip_clonalframe)
    """
    cp !{fasta} alignment-masked.fa
    pigz -c -n --best -p !{task.cpus} !{fasta} > !{params.prefix}.aligned.fa.gz
    """
    else
    """
    mkdir clonalframe

    iqtree -s !{fasta} -m !{params.m} -nt !{task.cpus} -fast -pre clonalframe/start-tree

    ClonalFrameML clonalframe/start-tree.treefile !{fasta} clonalframe/clonalframe \
        -emsim !{params.emsim} !{params.clonal_opts}

    maskrc-svg.py clonalframe/clonalframe --aln !{fasta} --symbol '-' \
        --out clonalframe/core_gene_alignment-masked.aln

    cp clonalframe/core_gene_alignment-masked.aln alignment-masked.fa
    pigz -c -n --best -p !{task.cpus} alignment-masked.fa > !{params.prefix}.aligned.fa.gz

    pigz -n --best -p !{task.cpus} clonalframe/core_gene_alignment-masked.aln
    """
}

process create_phylogeny {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "iqtree/*"
    publishDir outdir, mode: 'copy', overwrite: params.overwrite, pattern: "${params.prefix}.iqtree"

    input:
    file fasta from FINAL_TREE

    output:
    file 'iqtree/*'
    file "${params.prefix}.iqtree"

    shell:
    """
    mkdir iqtree
    iqtree -s !{fasta} -m !{params.m} -nt !{task.cpus} -pre iqtree/core-genome \
           -bb !{params.bb} -alrt !{params.alrt} -wbt -wbtl \
           -alninfo !{params.iqtree_opts}
    cp iqtree/core-genome.iqtree !{params.prefix}.iqtree
    """
}

process pairwise_snp_distance {
    publishDir outdir, mode: 'copy', overwrite: params.overwrite

    input:
    file fasta from SNP_DISTS

    output:
    file "${params.prefix}.distance.txt"

    shell:
    b = params.b ? "" : "-b"
    """
    snp-dists !{b} !{fasta} > !{params.prefix}.distance.txt
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
    error += is_positive_integer(params.i, 'i')
    error += is_positive_integer(params.cd, 'cd')
    error += is_positive_integer(params.g, 'g')
    error += is_positive_integer(params.emsim, 'emsim')
    error += is_positive_integer(params.bb, 'bb')
    error += is_positive_integer(params.alrt, 'alrt')
    error += is_positive_integer(params.prokka_coverage, 'prokka_coverage')

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


def check_gffs_exist(bactopia_path, sleep_time) {
    compressed = create_gff_channel(bactopia_path, false)
    log.info "\nIf this looks wrong, now's your chance to back out (CTRL+C 3 times)."
    log.info "Sleeping for ${sleep_time} seconds..."
    sleep(sleep_time * 1000)
    return compressed
}


def create_gff_channel(bactopia_path, is_process) {
    if (is_process) {
        return Channel.fromPath("${bactopia_path}/**/annotation/*.gff*")
    } else {
        gffs = Channel.fromPath("${bactopia_path}/**/annotation/*.gff*").toList()
        count = gffs.val.size()
        if (count > 0) {
            compressed = gffs.val[0].toString().endsWith("gz") ? true : false

            if (compressed) {
                log.info "Found ${count} compressed GFF files for analyis"
            } else {
                log.info "Found ${count} GFF files for analyis"
            }
            return compressed
        } else {
            log.error("Failed to find any GFF files. Please verify ${bactopia_path} contains Bactopia outputs.")
            exit 1
        }
    }
}


def print_help() {
    log.info"""
    Required Parameters:
        --bactopia STR          Directory containing Bactopia analysis results for all samples.

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

    RefSeq Assemblies Related Parameters:
        --species STR           The name of the species to download RefSeq assemblies for. This
                                    is a completely optional step and is meant to supplement
                                    your dataset with high-quality completed genomes.

        --prokka_evalue STR     Similarity e-value cut-off
                                    Default: ${params.prokka_evalue}

        --prokka_coverage INT   Minimum coverage on query protein
                                     Default: ${params.prokka_coverage}

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

    ClonalFrameML Related Parameters:
        --skip_clonalframe      Skip the ClonalFrameML and use the original core-genome
                                    alignment for the final tree.

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
    """.stripIndent()
    exit 0
}
