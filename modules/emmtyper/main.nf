process EMMTYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)
    path blastdb

    output:
    tuple val(meta), path("*.tsv")          , emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Conditionally add the database if it is provided by user
    if [ "$blastdb" == "" ]; then
        emmtyper \\
            $args \\
            $fasta_name \\
            > ${prefix}.tsv
    else
        # Make the blast database
        makeblastdb -in $blastdb -dbtype nucl

        emmtyper \\
            --blast_db $blastdb \\
            $args \\
            $fasta_name \\
            > ${prefix}.tsv

        # Remove the blast database
        rm $blastdb.*
    fi

    # If 'tmp' is not in $fasta_name, remove '.tmp' from the output files contents
    if [ $fasta_name != *.tmp* ]; then
        sed -i 's/.tmp\t/\t/g' ${prefix}.tsv
    fi


    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emmtyper: \$( echo \$(emmtyper --version 2>&1) | sed 's/^.*emmtyper v//' )
    END_VERSIONS
    """
}
