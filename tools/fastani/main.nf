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
samples = gather_sample_set(params.bactopia, params.exclude, params.include, params.sleep_time)
reference = null
if (params.reference) {
    if (file(params.reference).exists()) {
        samples << file(params.reference)
        reference = file(params.reference)
    } else {
        log.error("Could not open ${params.reference}, please verify existence. Unable to continue.")
        exit 1
    }
}


process collect_assemblies {    
    publishDir "${OUTDIR}/refseq", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "fasta/*.fna"
    publishDir "${OUTDIR}/refseq", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "accession-*.txt"

    input:
    file(fasta) from Channel.fromList(samples).collect()

    output:
    file 'fasta/*.fna' optional true
    file 'accession-*.txt' optional true
    file 'assemblies/*.fna' into ALL_FASTA
    file 'samples/*.fna.query' into EACH_SAMPLE

    shell:
    is_gzipped = null
    reference_name = null
    if (reference) {
        is_gzipped = params.reference.endsWith(".gz") ? true : false
        reference_name = reference.getSimpleName().replace(".gz", "")
        reference_file = reference.getName()
    }
    """
    mkdir assemblies samples
    if [ "!{params.species}" != "null" ]; then
        mkdir fasta
        if [ "!{params.limit}" != "null" ]; then
            ncbi-genome-download bacteria -l complete -o ./ -F fasta -p !{task.cpus} \
                                          --genus "!{params.species}" -r 50 --dry-run > accession-list.txt
            shuf accession-list.txt | head -n !{params.limit} > accession-subset.txt
            ncbi-genome-download bacteria -l complete -o ./ -F fasta -p !{task.cpus} \
                                          -A accession-subset.txt -r 50
        else
            ncbi-genome-download bacteria -l complete -o ./ -F fasta -p !{task.cpus} \
                                          --genus "!{params.species}" -r 50
        fi
        find -name "GCF*.fna.gz" | xargs -I {} mv {} fasta/
        rename 's/(GCF_\\d+).*/\$1.fna.gz/' fasta/*
        gunzip fasta/*
        cp fasta/*.fna assemblies/
        ls fasta/ | xargs -I {} cp -P fasta/{} samples/{}.query
    fi

    if [ "!{params.accession}" != "null" ]; then
        mkdir fasta
        ncbi-genome-download bacteria -l complete -o ./ -F fasta -p !{task.cpus} \
                                      -A !{params.accession} -r 50
        find -name "GCF*.fna.gz" | xargs -I {} mv {} fasta/
        rename 's/(GCF_\\d+).*/\$1.fna.gz/' fasta/*
        gunzip fasta/*
        cp fasta/*.fna assemblies/
        ls fasta/ | xargs -I {} cp -P fasta/{} samples/{}.query
    fi
    
    if [ "!{params.reference}" != "null" ]; then
        if [ "!{is_gzipped}" == "true" ]; then
            zcat !{reference_file} > assemblies/!{reference_name}.fna
            cp -P assemblies/!{reference_name}.fna samples/!{reference_name}.fna.query
        else
            cp -P !{reference_file} assemblies/!{reference_name}.fna
            cp -P !{reference_file} samples/!{reference_name}.fna.query
        fi
    fi

    find -name "*.fna.gz" | xargs -I {} gunzip {}

    if [ "!{params.pairwise}" == "true" ]; then
        cp -f -P *.fna assemblies/
    fi

    ls *.fna | xargs -I {} cp -f -P {} samples/{}.query
    """
}

process calculate_ani {
    tag "${query} - ${reference}"

    input:
    file(reference_fasta) from ALL_FASTA.flatten()
    each file(query_fasta) from EACH_SAMPLE

    output:
    file "${reference}.${query}.tsv" into RESULTS

    shell:
    reference = reference_fasta.getSimpleName()
    query = query_fasta.getSimpleName()
    """

    fastANI -q !{query_fasta} -r !{reference_fasta} --kmer !{params.kmer} --fragLen !{params.fragLen} \
            --minFraction !{params.minFraction} -o !{reference}.!{query}.tsv
    sed -i 's/.fna.query//' !{reference}.!{query}.tsv
    sed -i 's/.fna//' !{reference}.!{query}.tsv
    """
}

process merge_results {
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "fastani.tsv"
    publishDir OUTDIR, mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "references/*.tsv"

    input:
    file(results) from RESULTS.collect()

    output:
    file "fastani.tsv"
    file "references/*.tsv"

    shell:
    """
    mkdir references
    ls *.tsv | cut -f 1,1 -d "." | sort | uniq | \
        xargs -I {} sh -c 'cat {}.*.tsv | sort -k 2,2 -k3,3rn > references/{}.tsv' 
    cat *.tsv | sort -k 2,2 -k3,3rn > fastani.tsv
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
    IGNORE_LIST = ["container-path", "frag-len", "min-fraction"]
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
    error += is_positive_integer(params.kmer, 'kmer')
    error += is_positive_integer(params.fragLen, 'fragLen')
    error += is_positive_integer(params.limit, 'limit')

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

def build_assembly_tuple(sample, dir) {
    assembly = "${dir}/${sample}/assembly/${sample}.fna"
    if (file("${assembly}.gz").exists()) {
        // Compressed assemblies
        return file("${assembly}.gz")
    } else if (file(assembly).exists()) {
        return file(assembly)
    } else {
        log.error("Could not locate assembly for ${sample}, please verify existence. Unable to continue.")
        exit 1
    }
}

def gather_sample_set(bactopia_dir, exclude_list, include_list, sleep_time) {
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
    
    sample_list = []
    file(bactopia_dir).eachFile { item ->
        if( item.isDirectory() ) {
            sample = item.getName()
            if (!IGNORE_LIST.contains(sample)) {
                if (inclusions.contains(sample) || include_all) {
                    if (!exclusions.contains(sample)) {
                        if (is_sample_dir(sample, bactopia_dir)) {
                            sample_list << build_assembly_tuple(sample, bactopia_dir)
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
        This is a completely optional step and is meant to supplement your dataset with 
        high-quality completed genomes.

        --species STR           The name of the species to download RefSeq assemblies for.
        
        --accession STR         The Assembly accession (e.g. GCF*.*) download from RefSeq.

        --limit INT             Limit the number of RefSeq assemblies to download. If the the
                                    number of available genomes exceeds the given limit, a 
                                    random subset will be selected.
                                    Default: Download all available genomes        

    User Procided Reference:
        --reference STR         A reference genome to calculate 

    FastANI Related Parameters:
        --kmer INT              kmer size <= 16 
                                    Default: ${params.kmer}

        --fragLen INT           fragment length
                                    Default: ${params.fragLen}

        --minFraction FLOAT     Minimum fraction of genome that must be shared for trusting ANI.
                                    If reference and query genome size differ, smaller one among
                                    the two is considered.
                                    Default: ${params.minFraction}

        --pairwise              Pairwise ANI's will be calculated for all samples including any 
                                    downloaded RefSeq genomes or user provided genomes.


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
