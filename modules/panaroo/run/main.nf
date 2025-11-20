nextflow.preview.types = true

process PANAROO_RUN {
    tag "${prefix}"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, gff) : Tuple<Map, Set<Path>>

    stage:
    stageAs 'gff-tmp/*', gff

    output:
    supplemental = tuple(meta, file("panaroo/*"))
    aln          = tuple(meta, file("core-genome.aln.gz", optional: true))
    filtered_aln = tuple(meta, file("core-genome.filtered.aln.gz", optional: true))
    csv          = tuple(meta, file("panaroo/gene_presence_absence_roary.csv", optional: true))
    panaroo_csv  = tuple(meta, file("panaroo/gene_presence_absence.csv", optional: true))
    logs         = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin     = tuple(meta, file(".command.begin"))
    nf_err       = tuple(meta, file(".command.err"))
    nf_log       = tuple(meta, file(".command.log"))
    nf_out       = tuple(meta, file(".command.out"))
    nf_run       = tuple(meta, file(".command.run"))
    nf_sh        = tuple(meta, file(".command.sh"))
    nf_trace     = tuple(meta, file(".command.trace"))
    versions     = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "panaroo/logs"
    meta.process_name = task.ext.process_name
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gz" | xargs gunzip

    # Make FOFN of gff (Prokka) and gff3 (Bakta) files
    find gff/ -name "*.gff" -or -name "*.gff3" > gff-fofn.txt

    panaroo \\
        ${task.ext.args} \\
        -t ${task.cpus} \\
        -o supplemental \\
        -i gff-fofn.txt

    # Cleanup
    find . -name "*.fas" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.fa" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.fasta" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.aln" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.gml" | xargs -I {} -P ${task.cpus} -n 1 gzip {}

    if [[ -f "supplemental/core_gene_alignment.aln.gz" ]]; then
        mv supplemental/core_gene_alignment.aln.gz ./core-genome.aln.gz
    fi

    if [[ -f "supplemental/core_gene_alignment_filtered.aln.gz" ]]; then
        mv supplemental/core_gene_alignment_filtered.aln.gz ./core-genome.filtered.aln.gz
    fi

    if [[ -f "supplemental/gene_data.csv" ]]; then
        gzip supplemental/gene_data.csv
    fi

    mv supplemental/ panaroo/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //')
    END_VERSIONS
    """
}
