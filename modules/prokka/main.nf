process PROKKA {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("results/*.{fna,fna.gz}"), path("results/*.{faa,faa.gz}"), path("results/*.{gff,gff.gz}"), emit: annotations
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
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name

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
        mkdir tmp_prokka/
        TMPDIR=tmp_prokka/ bactopia-prokka \\
            $args \\
            --cpus $task.cpus \\
            --prefix $prefix \\
            ${compliant} \\
            ${locustag} \\
            $proteins_opt \\
            $prodigal_opt \\
            $fasta_name
        rm -rf tmp_prokka/
    else
        prokka \\
            $args \\
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

    # Cleanup intermediate files
    rm -rf ${fasta_name} blastdb/ results/*.pdb results/*.pjs results/*.pot results/*.ptf results/*.pto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makeblastdb: \$( echo \$(makeblastdb -version 2>&1) | sed 's/^.*makeblastdb: //;s/ .*\$//')
        prokka: \$( echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
