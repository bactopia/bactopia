process SHIGAPASS {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}_Flex_summary.tsv"), optional: true, emit: flex_tsv
    path "versions.yml", emit: versions
    path ".command.begin", emit: begin
    path ".command.err", emit: err
    path ".command.log", emit: log
    path ".command.out", emit: out
    path ".command.run", emit: run
    path ".command.sh", emit: sh
    path ".command.trace", emit: trace

    script:
    def args = task.ext.args ?: ''
    def SHIGAPASS_VERSION = "1.5.0"
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # ShigaPass does not accept compressed FASTA files
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Convert our genome path to a file with a path in it
    ls $fasta_name > ${fasta_name}_tmp.txt

    # Run ShigaPass
    ShigaPass.sh \\
        -l ${fasta_name}_tmp.txt \\
        $args \\
        -p "\$(dirname \$(which ShigaPass.sh))/../share/shigapass-${SHIGAPASS_VERSION}/db" \\
        -t $task.cpus \\
        -o ${prefix}

    # Remove the temporary file from above
    rm ${fasta_name}_tmp.txt ${fasta_name}

    # Convert to tab delimited and move to the pwd
    sed 's/;/\t/g' ${prefix}/ShigaPass_summary.csv > ${prefix}.tsv

    # Convert to tab delimited and move to the pwd
    [ ! -f ${prefix}/ShigaPass_Flex_summary.csv ] || sed 's/;/\t/g' ${prefix}/ShigaPass_Flex_summary.csv > ${prefix}_Flex_summary.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^.*ShigaPass version //' )
    END_VERSIONS
    """
}
