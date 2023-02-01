// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'bakta')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::bakta=1.5.1"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BAKTA_RUN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.6.0--pyhdfd78af_0' :
        'quay.io/biocontainers/bakta:1.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf
    path replicons

    output:
    tuple val(meta), path("results/${prefix}.{embl,embl.gz}")            , emit: embl
    tuple val(meta), path("results/${prefix}.{faa,faa.gz}")              , emit: faa
    tuple val(meta), path("results/${prefix}.{ffn,ffn.gz}")              , emit: ffn
    tuple val(meta), path("results/${prefix}.{fna,fna.gz}")              , emit: fna
    tuple val(meta), path("results/${prefix}.{gbff,gbff.gz}")            , emit: gbff
    tuple val(meta), path("results/${prefix}.{gff3,gff3.gz}")            , emit: gff
    tuple val(meta), path("results/${prefix}.hypotheticals.tsv")         , emit: hypotheticals_tsv
    tuple val(meta), path("results/${prefix}.hypotheticals.{faa,faa.gz}"), emit: hypotheticals_faa
    tuple val(meta), path("results/${prefix}.tsv")                       , emit: tsv
    tuple val(meta), path("results/${prefix}.txt")                       , emit: txt
    tuple val(meta), path("${prefix}-blastdb.tar.gz")                    , emit: blastdb
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions


    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    def replicons_opt = replicons ? "--replicons ${replicons[0]}" : ""
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$is_tarball" == "true" ]; then
        mkdir database
        tar -xzf $db -C database
        BAKTA_DB=\$(find database/ -name "bakta.db" | sed 's=bakta.db==')
    else
        BAKTA_DB=\$(find $db/ -name "bakta.db" | sed 's=bakta.db==')
    fi

    bakta \\
        --output results \\
        $options.args \\
        --threads $task.cpus \\
        --prefix ${prefix} \\
        --db \$BAKTA_DB \\
        $proteins_opt \\
        $prodigal_opt \\
        $replicons_opt \\
        $fasta

    # Make blastdb of contigs, genes, proteins
    mkdir blastdb
    cat ${prefix}/${prefix}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${prefix}" -out blastdb/${prefix}.fna
    cat ${prefix}/${prefix}.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for ${prefix}" -out blastdb/${prefix}.ffn
    cat ${prefix}/${prefix}.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for ${prefix}" -out blastdb/${prefix}.faa
    tar -czf ${prefix}-blastdb.tar.gz blastdb/

    if [[ "!{params.skip_compression}" == "false" ]]; then
        gzip --best results/${prefix}.embl
        gzip --best results/${prefix}.faa
        gzip --best results/${prefix}.ffn
        gzip --best results/${prefix}.fna
        gzip --best results/${prefix}.gbff
        gzip --best results/${prefix}.gff3
        gzip --best results/${prefix}.hypotheticals.faa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
        makeblastdb: $(echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    END_VERSIONS
    """
}

process BAKTA_MAIN_RUN {
    // This process is called in the main Bactopia pipeline
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.4.0--pyhdfd78af_1' :
        'quay.io/biocontainers/bakta:1.4.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(genome_size), path(fasta), path(total_contigs), path(proteins), path(prodigal_tf)
    path db
    path replicons

    output:
    tuple val(meta), path("results/${prefix}.{ffn,ffn.gz}"), path("results/${prefix}.{faa,faa.gz}"), emit: annotations
    tuple val(meta), path("results/${prefix}.{embl,embl.gz}")            , emit: embl
    tuple val(meta), path("results/${prefix}.{faa,faa.gz}")              , emit: faa
    tuple val(meta), path("results/${prefix}.{ffn,ffn.gz}")              , emit: ffn
    tuple val(meta), path("results/${prefix}.{fna,fna.gz}")              , emit: fna
    tuple val(meta), path("results/${prefix}.{gbff,gbff.gz}")            , emit: gbff
    tuple val(meta), path("results/${prefix}.{gff3,gff3.gz}")            , emit: gff
    tuple val(meta), path("results/${prefix}.hypotheticals.tsv")         , emit: hypotheticals_tsv
    tuple val(meta), path("results/${prefix}.hypotheticals.{faa,faa.gz}"), emit: hypotheticals_faa
    tuple val(meta), path("results/${prefix}.tsv")                       , emit: tsv
    tuple val(meta), path("results/${prefix}.txt")                       , emit: txt
    tuple val(meta), path("${prefix}-blastdb.tar.gz")                    , emit: blastdb
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions


    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def replicons_opt = replicons ? "--replicons ${replicons[0]}" : ""
    def prodigal_opt = ""
    if (prodigal_tf.getName() != 'EMPTY_TF' && !params.skip_prodigal_tf) {
        prodigal_opt = "--prodigaltf ${prodigal_tf}"
    }

    def proteins_opt = ""
    def genus = ""
    def species = ""
    if (proteins.getName() != 'EMPTY_PROTEINS') {
        proteins_opt = "--proteins ${proteins}"
        proteins_name = proteins.getName()
        if(proteins_name.contains("-")) {
            genus = "--genus ${proteins_name.split('-')[0].capitalize()}"
            species = "--species ${proteins_name.split('-')[1]}"
        } else {
            genus = "--genus ${proteins_name.capitalize()}"
            species = "--species spp."
        }
    }
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "$is_tarball" == "true" ]; then
        mkdir bakta
        tar -xzf $db -C bakta
        BAKTA_DB=\$(find bakta/ -name "bakta.db" | sed 's=bakta.db==')
    else
        BAKTA_DB=\$(find $db/ -name "bakta.db" | sed 's=bakta.db==')
    fi

    bakta \\
        --output results/ \\
        $options.args \\
        --threads $task.cpus \\
        --prefix ${prefix} \\
        --db \$BAKTA_DB \\
        $genus \\
        $species \\
        $proteins_opt \\
        $prodigal_opt \\
        $replicons_opt \\
        $fasta

    # Make blastdb of contigs, genes, proteins
    mkdir blastdb
    cat ${prefix}/${prefix}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${prefix}" -out blastdb/${prefix}.fna
    cat ${prefix}/${prefix}.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for ${prefix}" -out blastdb/${prefix}.ffn
    cat ${prefix}/${prefix}.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for ${prefix}" -out blastdb/${prefix}.faa
    tar -czf ${prefix}-blastdb.tar.gz blastdb/

    if [[ "${params.skip_compression}" == "false" ]]; then
        gzip --best results/${prefix}.embl
        gzip --best results/${prefix}.faa
        gzip --best results/${prefix}.ffn
        gzip --best results/${prefix}.fna
        gzip --best results/${prefix}.gbff
        gzip --best results/${prefix}.gff3
        gzip --best results/${prefix}.hypotheticals.faa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
        makeblastdb: $(echo $(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*$//')
    END_VERSIONS
    """
}
