#! /usr/bin/env nextflow
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version
OUTDIR = "${params.outdir}/bactopia-tools/${PROGRAM_NAME}/${params.prefix}"
OVERWRITE = workflow.resume || params.force ? true : false
DUMMY_NAME = "DUMMY_FILE"

// Validate parameters
if (params.version) print_version();
log.info "bactopia tools ${PROGRAM_NAME} - ${VERSION}"
if (params.help || workflow.commandLine.trim().endsWith(workflow.scriptName)) print_help();
check_input_params()
samples = gather_sample_set(params.bactopia, params.exclude, params.include, params.sleep_time, params.only_completed)
accessions = [tuple(null, null)]
if (params.accessions) {
    if (file(params.accessions).exists()) {
        accessions = [tuple(true, file(params.accessions))]
    } else {
        log.error("Could not open ${params.accessions}, please verify existence. Unable to continue.")
        exit 1
    }
}

process download_references {
    publishDir "${OUTDIR}/refseq", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "fasta/*.fna"
    publishDir "${OUTDIR}/refseq", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "accession-*.txt"

    input:
    set val(has_accessions), file(accession_list) from accessions

    output:
    file("fasta/*.fna") into ANNOTATE
    file 'accession-*.txt' optional true

    shell:
    opt = false
    if (params.species) {
        opt = "-g '${params.species}'"
    } else if (params.accession) {
        opt = "-A '${params.accession}'"
    }
    
    """
    mkdir fasta
    if [ "!{opt}" == "false" ]; then
        touch fasta/!{DUMMY_NAME}.fna
    else
        if [ "!{params.species}" != "null" ]; then
            if [ "!{params.limit}" != "null" ]; then
                ncbi-genome-download bacteria -l complete -o ./ -F fasta \
                                              !{opt} -r 50 --dry-run > accession-list.txt
                shuf accession-list.txt | head -n !{params.limit} | cut -f 1,1 > accession-subset.txt
                ncbi-genome-download bacteria -l complete -o ./ -F fasta \
                                              -A accession-subset.txt -r 50
            else
                ncbi-genome-download bacteria -l complete -o ./ -F fasta -p !{task.cpus} !{opt} -r 50
            fi
        elif [ "!{has_accessions}" == "true" ]; then
            ncbi-genome-download bacteria -l complete -o ./ -F fasta -A accession_list -r 50
        else
            ncbi-genome-download bacteria -l complete -o ./ -F fasta !{opt} -r 50
        fi
        find -name "GCF*.fna.gz" | xargs -I {} mv {} fasta/
        rename 's/(GCF_\\d+).*/\$1.fna.gz/' fasta/*
        gunzip fasta/*
    fi
    """
}

process annotate_references {
    publishDir "${OUTDIR}/refseq", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "gff/${name}.gff"
    tag "${name}"

    input:
    file(fasta) from ANNOTATE.flatten()

    output:
    file "gff/${name}.gff" into REFERENCE_GFF

    shell:
    name = fasta.getSimpleName()
    """
    if [ "!{name}" == "!{DUMMY_NAME}" ]; then
        mkdir gff
        touch gff/!{DUMMY_NAME}.gff
    else
        prokka !{fasta} --cpus !{task.cpus} --evalue '!{params.prokka_evalue}' \
                        --outdir gff --prefix !{name} --coverage !{params.prokka_coverage}
    fi
    """

}

process build_pangenome {
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "roary/*"

    input:
    file(sample_gff) from Channel.fromList(samples).collect().ifEmpty { "${DUMMY_NAME}.gff" }
    file(reference_gff) from REFERENCE_GFF.collect()

    output:
    file 'roary/*'
    file 'alignment.fa' into RECOMBINATION

    shell:
    n = params.n ? "-n" : ""
    s = params.s ? "-s" : ""
    ap = params.ap ? "-ap" : ""
    """
    rm -rf !{DUMMY_NAME}.gff
    find -name "*.gff.gz" | xargs -I {} gunzip {}
    roary -f roary -e !{n} -v -p !{task.cpus} !{s} !{ap} -g !{params.g} \
          -i !{params.i} -cd !{params.cd} -iv !{params.iv} -r *.gff
    cp roary/core_gene_alignment.aln alignment.fa
    pigz -n --best -p !{task.cpus} roary/core_gene_alignment.aln
    """
}

process identify_recombination {
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "clonalframe/*"
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "${params.prefix}.aligned.fa.gz"

    input:
    file fasta from RECOMBINATION

    output:
    file 'clonalframe/*' optional true
    file "${params.prefix}.aligned.fa.gz"
    file 'alignment-masked.fa' into FINAL_TREE, SNP_DISTS

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
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "iqtree/*"
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "${params.prefix}.iqtree"

    input:
    file fasta from FINAL_TREE

    output:
    file 'iqtree/*'
    file "${params.prefix}.iqtree"

    when:
    params.skip_phylogeny == false

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
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE

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
    error += is_positive_integer(params.i, 'i')
    error += is_positive_integer(params.cd, 'cd')
    error += is_positive_integer(params.g, 'g')
    error += is_positive_integer(params.emsim, 'emsim')
    error += is_positive_integer(params.bb, 'bb')
    error += is_positive_integer(params.alrt, 'alrt')
    error += is_positive_integer(params.prokka_coverage, 'prokka_coverage')

    if (params.limit) {
        error += is_positive_integer(params.limit, 'limit')
    }

    if (params.only_completed && !params.species && !params.accession) {
        log.error("'--only_completed' requires that '--species' or '--accession' is used also.")
        error += 1
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


def is_sample_dir(sample, dir){
    return file("${dir}/${sample}/${sample}-genome-size.txt").exists()
}

def build_gff_tuple(sample, dir) {
    gff = "${dir}/${sample}/annotation/${sample}.gff"
    if (file("${gff}.gz").exists()) {
        // Compressed assemblies
        return file("${gff}.gz")
    } else if (file(gff).exists()) {
        return file(gff)
    } else {
        log.error("Could not locate GFF for ${sample}, please verify existence. Unable to continue.")
        exit 1
    }
}

def gather_sample_set(bactopia_dir, exclude_list, include_list, sleep_time, only_completed) {
    sample_list = []
    if (only_completed) {
        log.info "--only_completed given, so only the completed genomes will be used for pan-genome analysis."
    } else {
        include_all = true
        inclusions = []
        exclusions = []
        IGNORE_LIST = ['.nextflow', 'bactopia-info', 'bactopia-tools', 'work',]
        if (include_list) {
            new File(include_list).eachLine { line -> 
                inclusions << line.trim().split('\t')[0]
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
        
        file(bactopia_dir).eachFile { item ->
            if( item.isDirectory() ) {
                sample = item.getName()
                if (!IGNORE_LIST.contains(sample)) {
                    if (inclusions.contains(sample) || include_all) {
                        if (!exclusions.contains(sample)) {
                            if (is_sample_dir(sample, bactopia_dir)) {
                                sample_list << build_gff_tuple(sample, bactopia_dir)
                            } else {
                                log.info "${sample} is missing genome size estimate file"
                            }
                        }
                    }
                }
            }
        }

        if (sample_list.size == 0) {
            log.error "Did not find any samples in ${bactopia_dir} to process, please check and try again."
            exit 1
        } else {
            log.info "Found ${sample_list.size} samples to process"
        }
    }

    log.info "\nIf this looks wrong, now's your chance to back out (CTRL+C 3 times)."
    log.info "Sleeping for ${sleep_time} seconds..."
    sleep(sleep_time * 1000)
    return sample_list
}

def print_help() {
    log.info"""
    Required Parameters:
        --bactopia STR          Directory containing Bactopia analysis results for all samples.

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

    RefSeq Assemblies Related Parameters:
        --species STR           The name of the species to download RefSeq assemblies for. This
                                    is a completely optional step and is meant to supplement
                                    your dataset with high-quality completed genomes.

        --accession STR         A NCBI Assembly database RefSeq accession to be downloaded and included
                                    in the pan-genome analysis.

        --accessions STR        A file with Assembly accessions (e.g. GCF*.*) to download from RefSeq.

        --limit INT             Limit the number of RefSeq assemblies to download. If the the
                                    number of available genomes exceeds the given limit, a 
                                    random subset will be selected.
                                    Default: Download all available genomes
        
        --only_completed        Pan-genome will be created using only the completed RefSeq genomes.    

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
