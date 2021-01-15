#! /usr/bin/env nextflow
PROGRAM_NAME = workflow.manifest.name
VERSION = workflow.manifest.version
OUTDIR = "${params.outdir}/bactopia-tools/${PROGRAM_NAME}/${params.prefix}"
DOWNLOAD_EGGNOG = false
OVERWRITE = workflow.resume || params.force ? true : false

// Validate parameters
if (params.version) print_version();
log.info "bactopia tools ${PROGRAM_NAME} - ${VERSION}"
if (params.help || workflow.commandLine.trim().endsWith(workflow.scriptName)) print_help();
check_input_params()
samples = gather_sample_set(params.bactopia, params.exclude, params.include, params.sleep_time)
EMAPPER_PARAMS = build_emapper_params()
log.info "EMAPPER_PARAMS -> ${EMAPPER_PARAMS}"

process download_eggnogdb {
    publishDir "${params.eggnog}", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "eggnog.db"
    publishDir "${params.eggnog}", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "eggnog_proteins.dmnd"

    input:
    val eggnog_path from params.eggnog

    output:
    file("eggnog.db") optional true
    file("eggnog_proteins.dmnd") optional true

    when:
    DOWNLOAD_EGGNOG == true

    shell:
    """
    download_eggnog_data.py -y --data_dir ./ && touch completed.txt
    """
}

process eggnog_mapper {
    publishDir "${OUTDIR}", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "${sample}.emapper.annotations"
    publishDir "${OUTDIR}", mode: "${params.publish_mode}", overwrite: OVERWRITE, pattern: "${sample}.emapper.seed_orthologs"
    tag "${sample}"

    input:
    file(eggnogdb) from Channel.fromPath( "${params.eggnog}/*" ).collect()
    each file(fasta) from Channel.fromList(samples)
    
    output:
    file("${sample}.emapper.annotations")
    file("${sample}.emapper.seed_orthologs")

    shell:
    sample = fasta.getSimpleName()
    is_gzipped = fasta.getName().endsWith('gz') ? true : false
    """
    if [ "!{is_gzipped}" == "true" ]; then
        zcat !{fasta} > !{sample}.faa
    fi
    emapper.py -i !{sample}.faa \
        --output_dir ./ \
        --output !{sample} \
        !{EMAPPER_PARAMS} -m diamond \
        --cpu !{task.cpus} \
        --data_dir ./ 
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

def build_emapper_params() {
    // String together the emapper.py params
    emapper_params = [
        params.keep_mapping_files ? '--keep_mapping_files' : '',
        params.no_annot ? '--no_annot' : '',
        params.no_file_comments ? '--no_file_comments' : '',
        params.no_refine ? '--no_refine' : '',
        params.no_search ? '--no_refine' : '',
        params.tax_scope ? "--tax_scope ${params.tax_scope}" : '',
        params.go_evidence ? "--go_evidence ${params.go_evidence}" : '',
        params.hmm_maxhits ? "--hmm_maxhits ${params.hmm_maxhits}" : '',
        params.hmm_evalue ? "--hmm_evalue ${params.hmm_evalue}" : '',
        params.hmm_score ? "--hmm_score ${params.hmm_score}" : '',
        params.hmm_maxseqlen ? "--hmm_maxseqlen ${params.hmm_maxseqlen}" : '',
        params.matrix ? "--matrix ${params.matrix}" : '',
        params.gapopen ? "--gapopen ${params.gapopen}" : '',
        params.gapextend ? "--matrix ${params.gapextend}" : '',
        params.query_cover ? "--query-cover ${params.query_cover}" : '',
        params.subject_cover ? "--subject-cover ${params.subject_cover}" : '',
        params.seed_ortholog_evalue ? "--seed_ortholog_evalue ${params.seed_ortholog_evalue}" : '',
        params.seed_ortholog_score ? "--seed_ortholog_score ${params.seed_ortholog_score}" : '',
        params.target_taxa ? "--target_taxa ${params.target_taxa}" : '',
        params.predict_output_format ? "--predict_output_format ${params.predict_output_format}" : '',
    ]
    return (emapper_params.join(' ').replaceAll("\\s{2,}", " ").trim())
}

def check_input_params() {
    // Check for unexpected paramaters
    error = check_unknown_params()

    if (params.bactopia) {
        error += file_exists(params.bactopia, '--bactopia')
    } 

    if (params.eggnog) {
        eggnog_db = "${params.eggnog}/eggnog.db"
        diamond_db = "${params.eggnog}/eggnog_proteins.dmnd"
        if (file(eggnog_db).exists() && params.download_eggnog) {
            log.info """
                Found eggNOG database at ${params.eggnog}, but '--download_eggnog'
                also given. Existing eggNOG database will be over written with
                latest build. If this is not wanted, please stop now.
            """.stripIndent()
            DOWNLOAD_EGGNOG = true
        } else if (!file(eggnog_db).exists() && params.download_eggnog) {
            log.info("Latest eggNOG database will be downloaded to ${params.eggnog}")
            DOWNLOAD_EGGNOG = true
        } else if (!file(eggnog_db).exists()) {
            log.error """
                Please check that a eggNOG database exists at ${params.eggnog}. 
                Otherwise use '--download_eggnog' to download the latest eggNOG database
            """.stripIndent()
            error += 1
        }
    }

    if (params.include) {
        error += file_exists(params.include, '--include')
    } 
    
    if (params.exclude) {
        error += file_exists(params.exclude, '--exclude')
    } 
    
    if (!params.bactopia || !params.eggnog) {
        log.error """
        The required '--bactopia' and '--eggnog' parameters are missing, please check and try again.

        Required Parameters:
            --bactopia STR          Directory containing Bactopia analysis results for all samples.

            --eggnog STR            Directory containing the the eggNOG database files:
                                        eggnog.db and eggnog_proteins.dmnd.  If the database is not 
                                        found, you must use '--download_eggnog'.
        """.stripIndent()
        error += 1
    }

    error += is_positive_integer(params.cpus, 'cpus')
    error += is_positive_integer(params.max_time, 'max_time')
    error += is_positive_integer(params.max_memory, 'max_memory')
    error += is_positive_integer(params.sleep_time, 'sleep_time')

    error += is_positive_integer(params.hmm_maxhits, 'hmm_maxhits')
    error += is_positive_integer(params.hmm_evalue, 'hmm_evalue')
    error += is_positive_integer(params.hmm_score, 'hmm_score')
    error += is_positive_integer(params.hmm_maxseqlen, 'hmm_maxseqlen')
    error += is_positive_integer(params.Z, 'Z')
    error += is_positive_integer(params.query_cover, 'query_cover')
    error += is_positive_integer(params.subject_cover, 'subject_cover')
    error += is_positive_integer(params.seed_ortholog_evalue, 'seed_ortholog_evalue')
    error += is_positive_integer(params.seed_ortholog_score, 'seed_ortholog_score')

    if (params.gapopen) {
        error += is_positive_integer(params.gapopen, 'gapopen')
    }

    if (params.gapextend) {
        error += is_positive_integer(params.gapextend, 'gapextend')
    }

    if (params.hmm_qcov) {
        error += is_positive_integer(params.hmm_qcov, 'hmm_qcov') 
    }

    // Check for valid choices
    if (params.target_orthologs) {
        ALLOWED_TARGETS = ['one2one', 'many2one', 'one2many', 'many2many', 'all']
        if (!ALLOWED_TARGETS.contains(params.target_orthologs)) {
            log.error("'${params.target_orthologs}' is not a valid coice for --target_orthologs. Choices are: ${ALLOWED_TARGETS}")
            error += 1
        }
    }

    if (params.go_evidence) {
        ALLOWED_GO = ['experimental', 'non-electronic']
        if (!ALLOWED_GO.contains(params.go_evidence)) {
            log.error("'${params.go_evidence}' is not a valid coice for --go_evidence. Choices are: ${ALLOWED_GO}")
            error += 1
        }
    }

    if (params.matrix) {
        ALLOWED_MATRIX = ['BLOSUM62', 'BLOSUM90', 'BLOSUM80', 'BLOSUM50', 'BLOSUM45', 'PAM250', 'PAM70', 'PAM30']
        if (!ALLOWED_MATRIX.contains(params.matrix)) {
            log.error("'${params.matrix}' is not a valid coice for --matrix. Choices are: ${ALLOWED_MATRIX}")
            error += 1
        }
    }

    if (params.predict_output_format) {
        ALLOWED_FORMAT = ['per_query', 'per_species']
        if (!ALLOWED_FORMAT.contains(params.predict_output_format)) {
            log.error("'${params.predict_output_format}' is not a valid coice for --predict_output_format. Choices are: ${ALLOWED_FORMAT}")
            error += 1
        }
    }

    if (!['dockerhub', 'github', 'quay'].contains(params.registry)) {
            log.error "Invalid registry (--registry ${params.registry}), must be 'dockerhub', " +
                      "'github' or 'quay'. Please correct to continue."
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
    if (value) {
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
    }
    return error
}


def is_sample_dir(sample, dir){
    return file("${dir}/${sample}/${sample}-genome-size.txt").exists()
}

def build_protein_tuple(sample, dir) {
    proteins = "${dir}/${sample}/annotation/${sample}.faa"
    if (file("${proteins}.gz").exists()) {
        // Compressed proteins
        tuple(file("${proteins}.gz"))
    } else if (file(proteins).exists()) {
        tuple(file(proteins))
    } else {
        log.error("Could not locate proteins for ${sample}, please verify existence. Unable to continue.")
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
                            sample_list << build_protein_tuple(sample, bactopia_dir)
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

        --eggnog STR            Directory containing the the eggNOG database files:
                                    eggnog.db and eggnog_proteins.dmnd.  If the database is not 
                                    found, you must use '--download_eggnog'.
                                    WARNING: eggNOG databases stored on NFS will see a significant
                                             increase in runtimes. If possible, SSD or Ramdisk 
                                             should be used.

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

    eggNOG-mapper Parameters:
    *Note: Unless specified the eggNOG-mapper defaults are used.*

        --download_eggnog       Download the latest eggNOG database, even if it exists.

    eggNOG Annotation Parameters:
        --tax_scope STR         Fix the taxonomic scope used for annotation, so only orthologs
                                    from a particular clade are used for functional transfer. 
        
        --target_orthologs STR  Defines what type of orthologs should be used for functional 
                                    transfer.
                                    Choices are: one2one, many2one, one2many, many2many, all

        --go_evidence STR       Defines what type of GO terms should be used for annotation. Choices are:
                                    'experimental': Use only terms inferred from experimental evidence
                                    'non-electronic': Use only non- electronically curated terms

    eggNOG HMM Search Parameters:
        --hmm_maxhits INT       Max number of hits to report.

        --hmm_evalue FLOAT      E-value threshold.

        --hmm_score INT         Bit score threshold.

        --hmm_maxseqlen INT     Ignore query sequences larger than `maxseqlen`.

        --hmm_qcov FLOAT        Min query coverage (from 0 to 1).

        --Z INT                 Fixed database size used in phmmer/hmmscan (allows comparing e-values 
                                    among databases).

    eggNOG DIAMOND Search Parameters:
        --use_diamond           Use DIAMOND instead of HMMER.

        --matrix STR            Scoring matrix. Choices are: BLOSUM62, BLOSUM90, BLOSUM80, BLOSUM50,
                                                             BLOSUM45, PAM250, PAM70, PAM30
        
        --gapopen INT           Gap open penalty

        --gapextend INT         Gap extend penalty

        --query_cover FLOAT     Report only alignments above the given percentage of query cover.

        --subject_cover FLOAT   Report only alignments above the given percentage of subject cover.

    eggNOG Seed Ortholog Search Parameters:
        --seed_ortholog_evalue FLOAT
                                Min E-value expected when searching for seed eggNOG ortholog. 
                                    Applies to phmmer/diamond searches. Queries not having a 
                                    --significant seed orthologs will not be annotated. 

        --seed_ortholog_score INT
                                Min bit score expected when searching for seed eggNOG ortholog. 
                                    Applies to phmmer/diamond searches. Queries not having a 
                                    --significant seed orthologs will not be annotated. 

    eggNOG Output Parameters:
        --keep_mapping_files    Do not delete temporary mapping files used for annotation (i.e. 
                                    HMMER and DIAMOND search outputs)

        --no_annot              Skip functional annotation, reporting only hits

        --no_file_comments      No header lines nor stats are included in the output files

        --no_refine             Skip hit refinement, reporting only HMM hits.

        --no_search             Skip HMM search mapping. Use existing hits file

    eggNOG Predict Orthologs Parameters:
        --target_taxa STR       Taxa that will be searched for orthologs

        --predict_output_format STR
                                Choose the output format among: per_query, per_species. 

    Nextflow Related Parameters:
        --condadir DIR          Directory to Nextflow should use for Conda environments
                                    Default: Bactopia's Nextflow directory


        --registry STR          Docker registry to pull containers from. 
                                    Available options: dockerhub, quay, or github
                                    Default: dockerhub

        --singularity_cache STR Directory where remote Singularity images are stored. If using a cluster, it must
                                    be accessible from all compute nodes.
                                    Default: NXF_SINGULARITY_CACHEDIR evironment variable, otherwise ${params.singularity_cache}


        --cleanup_workdir       After Bactopia is successfully executed, the work directory will be deleted.
                                    Warning: by doing this you lose the ability to resume workflows.
                                    
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
