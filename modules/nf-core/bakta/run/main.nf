// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'bakta')
options.btype = "main"
conda_tools   = "bioconda::bakta=1.9.4"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BAKTA_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.9.4--pyhdfd78af_0' :
        'quay.io/biocontainers/bakta:1.9.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf
    path replicons

    output:
    tuple val(meta), path("results/${prefix}.{fna,fna.gz}"), path("results/${prefix}.{faa,faa.gz}"), path("results/${prefix}.{gff3,gff3.gz}"), emit: annotations
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
        BAKTA_DB=\$(find database/ -name "bakta.db" | sed 's=bakta\\.db==')
    else
        BAKTA_DB=\$(find $db/ -name "bakta.db" | sed 's=bakta\\.db==')
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
    cat results/${prefix}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${prefix}" -out blastdb/${prefix}.fna
    cat results/${prefix}.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for ${prefix}" -out blastdb/${prefix}.ffn
    cat results/${prefix}.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for ${prefix}" -out blastdb/${prefix}.faa
    tar -cvf - blastdb/ | gzip -c > ${prefix}-blastdb.tar.gz

    if [[ "${params.skip_compression}" == "false" ]]; then
        gzip --best results/${prefix}.embl
        gzip --best results/${prefix}.faa
        gzip --best results/${prefix}.ffn
        gzip --best results/${prefix}.fna
        gzip --best results/${prefix}.gbff
        gzip --best results/${prefix}.gff3
        gzip --best results/${prefix}.hypotheticals.faa
    fi

    # Clean up
    if [ "$is_tarball" == "true" ]; then
        rm -rf database
    fi
    rm -rf blastdb/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
        makeblastdb: \$(echo \$(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*\$//')
    END_VERSIONS
    """
}
