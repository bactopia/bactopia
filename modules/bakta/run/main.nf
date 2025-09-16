process BAKTA_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf
    path replicons

    output:
    tuple val(meta), path("bakta/${prefix}.{fna,fna.gz}"), path("bakta/${prefix}.{faa,faa.gz}"), path("bakta/${prefix}.{gff3,gff3.gz}"), emit: annotations
    tuple val(meta), path("bakta/${prefix}.{embl,embl.gz}")            , emit: embl
    tuple val(meta), path("bakta/${prefix}.{faa,faa.gz}")              , emit: faa
    tuple val(meta), path("bakta/${prefix}.{ffn,ffn.gz}")              , emit: ffn
    tuple val(meta), path("bakta/${prefix}.{fna,fna.gz}")              , emit: fna
    tuple val(meta), path("bakta/${prefix}.{gbff,gbff.gz}")            , emit: gbff
    tuple val(meta), path("bakta/${prefix}.{gff3,gff3.gz}")            , emit: gff
    tuple val(meta), path("bakta/${prefix}.hypotheticals.tsv")         , emit: hypotheticals_tsv
    tuple val(meta), path("bakta/${prefix}.hypotheticals.{faa,faa.gz}"), emit: hypotheticals_faa
    tuple val(meta), path("bakta/${prefix}.tsv")                       , emit: tsv
    tuple val(meta), path("bakta/${prefix}.txt")                       , emit: txt
    tuple val(meta), path("bakta/${prefix}-blastdb.tar.gz")            , emit: blastdb
    tuple val(meta), path("*.{log,err}")                               , emit: logs, optional: true
    tuple val(meta), path(".command.begin")                            , emit: nf_begin
    tuple val(meta), path(".command.err")                              , emit: nf_err
    tuple val(meta), path(".command.log")                              , emit: nf_log
    tuple val(meta), path(".command.out")                              , emit: nf_out
    tuple val(meta), path(".command.run")                              , emit: nf_run
    tuple val(meta), path(".command.sh")                               , emit: nf_sh
    tuple val(meta), path(".command.trace")                            , emit: nf_trace
    tuple val(meta), path("versions.yml")                              , emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/main/annotator/bakta"
    meta.logs_dir = "${meta.id}/main/annotator/bakta/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
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
        --output bakta \\
        $task.ext.args \\
        --threads $task.cpus \\
        --prefix ${prefix} \\
        --db \$BAKTA_DB \\
        $proteins_opt \\
        $prodigal_opt \\
        $replicons_opt \\
        $fasta

    # Make blastdb of contigs, genes, proteins
    mkdir blastdb
    cat bakta/${prefix}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${prefix}" -out blastdb/${prefix}.fna
    cat bakta/${prefix}.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for ${prefix}" -out blastdb/${prefix}.ffn
    cat bakta/${prefix}.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for ${prefix}" -out blastdb/${prefix}.faa
    tar -cvf - blastdb/ | gzip -c > bakta/${prefix}-blastdb.tar.gz

    if [[ "${task.ext.skip_compression}" == "false" ]]; then
        gzip --best bakta/${prefix}.embl
        gzip --best bakta/${prefix}.faa
        gzip --best bakta/${prefix}.ffn
        gzip --best bakta/${prefix}.fna
        gzip --best bakta/${prefix}.gbff
        gzip --best bakta/${prefix}.gff3
        gzip --best bakta/${prefix}.hypotheticals.faa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
