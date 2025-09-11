process CLERMONTYPING {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv") , emit: tsv
    tuple val(meta), path("supplemental/*"), emit: results
    tuple val(meta), path("*.{log,err}"   ), emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    clermonTyping.sh \\
        --fasta $fasta_name \\
        --name supplemental \\
        $task.ext.args

    # Remove temporary files and rename outputs
    rm supplemental/${fasta_name} ${fasta_name}
    rm -rf supplemental/db
    rm -rf supplemental/supplemental.R
    mv supplemental/${fasta_name}.xml supplemental/${prefix}.blast.xml
    mv supplemental/${fasta_name}_mash_screen.tab supplemental/${prefix}.mash.tsv
    mv supplemental/supplemental.html supplemental/${prefix}.html

    # add column names to phylogroups file
    echo "sample<TAB>detected_genes<TAB>quadruplex_genes<TAB>quadruplex_alleles<TAB>phylogroup<TAB>mash_group" | sed 's/<TAB>/\t/g' > ${prefix}.tsv
    cat supplemental/supplemental_phylogroups.txt >> ${prefix}.tsv
    rm supplemental/supplemental_phylogroups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clermontyping: \$(echo \$(clermonTyping.sh -v 2>&1) | sed 's/^.* version : //;s/ .*\$//')
    END_VERSIONS
    """
}
