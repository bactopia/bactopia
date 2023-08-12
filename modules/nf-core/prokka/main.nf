// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'prokka')
options.btype = options.btype ?: "main"
conda_tools   = "bioconda::prokka=1.14.6 bioconda::perl-bioperl=1.7.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PROKKA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0' :
        'quay.io/biocontainers/prokka:1.14.6--pl526_0' }"

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("results/*.{ffn,ffn.gz}"), path("results/*.{faa,faa.gz}"), emit: annotations
    tuple val(meta), path("results/*.{gff,gff.gz}"), emit: gff
    tuple val(meta), path("results/*.{gbk,gbk.gz}"), emit: gbk
    tuple val(meta), path("results/*.{fna,fna.gz}"), emit: fna
    tuple val(meta), path("results/*.{faa,faa.gz}"), emit: faa
    tuple val(meta), path("results/*.{ffn,ffn.gz}"), emit: ffn
    tuple val(meta), path("results/*.{sqn,sqn.gz}"), emit: sqn
    tuple val(meta), path("results/*.{fsa,fsa.gz}"), emit: fsa
    tuple val(meta), path("results/*.{tbl,tbl.gz}"), emit: tbl
    tuple val(meta), path("results/*.txt")         , emit: txt
    tuple val(meta), path("results/*.tsv")         , emit: tsv
    tuple val(meta), path("results/${prefix}-blastdb.tar.gz"), emit: blastdb
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")

    // Contig ID must <= 37 characters
    def compliant = params.compliant ? "--compliant" : ""
    def locustag = "--locustag ${meta.id}"
    if ("gnl|${params.centre}|${meta.id}_100000".length() > 37) {
        locustag = ""
        compliant = "--compliant"
    }
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    if [ "$params.prokka_debug" == "true" ]; then
        export PROKKA_DBDIR=\$(echo "\$(which prokka | sed "s=/prokka==")/../db")
        env
        bactopia-prokka \\
            $options.args \\
            --cpus $task.cpus \\
            --prefix $prefix \\
            ${compliant} \\
            ${locustag} \\
            $proteins_opt \\
            $prodigal_opt \\
            $fasta_name
    else
        prokka \\
            $options.args \\
            --cpus $task.cpus \\
            --prefix $prefix \\
            ${compliant} \\
            ${locustag} \\
            $proteins_opt \\
            $prodigal_opt \\
            $fasta_name
    fi

    # Make blastdb of contigs, genes, proteins
    mkdir blastdb
    cat ${prefix}/${prefix}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${prefix}" -out blastdb/${prefix}.fna
    cat ${prefix}/${prefix}.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for ${prefix}" -out blastdb/${prefix}.ffn
    cat ${prefix}/${prefix}.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for ${prefix}" -out blastdb/${prefix}.faa
    tar -cvf - blastdb/ | gzip -c > ${prefix}/${prefix}-blastdb.tar.gz

    if [[ "${params.skip_compression}" == "false" ]]; then
        gzip ${prefix}/*.gff
        gzip ${prefix}/*.gbk
        gzip ${prefix}/*.fna
        gzip ${prefix}/*.faa
        gzip ${prefix}/*.ffn
        gzip ${prefix}/*.sqn
        gzip ${prefix}/*.fsa
        gzip ${prefix}/*.tbl
    fi

    mv ${prefix}/${prefix}.err ./
    mv ${prefix}/${prefix}.log ./
    mv ${prefix}/ results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makeblastdb: \$( echo \$(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*\$//')
        prokka: \$( echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
