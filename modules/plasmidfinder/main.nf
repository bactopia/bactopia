nextflow.preview.types = true

process PLASMIDFINDER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    output:
    json        = tuple(meta, files("*.json"))
    txt         = tuple(meta, files("*.txt"))
    tsv         = tuple(meta, file("${prefix}.tsv"))
    genome_seq  = tuple(meta, files("*-hit_in_genome_seq.fsa.gz"))
    plasmid_seq = tuple(meta, files("*-plasmid_seqs.fsa.gz"))
    logs        = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs     = tuple(meta, files(".command.*"))
    versions    = tuple(meta, file("versions.yml"))

    script:
    def VERSION = '2.1.6'
    // Version information not provided by tool on CLI
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    plasmidfinder.py \\
        ${task.ext.args} \\
        -i ${fasta_name} \\
        -o ./ \\
        -x

    # Rename hard-coded outputs with prefix to avoid name collisions
    mv data.json ${prefix}.json
    mv results.txt ${prefix}.txt
    mv Hit_in_genome_seq.fsa ${prefix}-hit_in_genome_seq.fsa
    mv Plasmid_seqs.fsa ${prefix}-plasmid_seqs.fsa

    # Add sample name to TSV results
    head -n 1 results_tab.tsv | sed "s/^/Sample\t/" > ${prefix}.tsv
    tail -n +2 results_tab.tsv | sed "s/^/${prefix}\t/" >> ${prefix}.tsv

    # Cleanup
    gzip *.fsa
    rm -rf ${fasta_name} results_tab.tsv tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: ${VERSION}
    END_VERSIONS
    """
}
