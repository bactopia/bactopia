/**
 * A fast and scalable tool for bacterial pangenome analysis.
 *
 * This process executes panaroo_run to perform analysis
 *
 * @status stable
 * @keywords gff, pan-genome, alignment
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output, compression, conditional-logic
 * @citation panaroo_run
 *
 * @input tuple(meta, gff)
 * - `meta`: Groovy Map containing sample information
 * - `gff`: A set of GFF3 formatted files
 *
 * @output supplemental Supplemental
 * @output aln          Core-genome alignment produced by Panaroo (Optional)
 * @output filtered_aln Filtered Aln
 * @output csv          Gene presence absence in Roary format (Optional)
 * @output panaroo_csv  Gene presence absence in Panaroo format (Optional)
 * @output logs         Optional tool execution logs
 * @output nf_logs      Nextflow execution logs
 * @output versions     Software version information (YAML format)
 */
nextflow.preview.types = true

process PANAROO_RUN {
    tag "${prefix}"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, gff) : Tuple<Map, List<Path>>

    stage:
    stageAs 'gff-tmp/*', gff

    output:
    supplemental = tuple(meta, files("panaroo/*"))
    aln          = tuple(meta, file("core-genome.aln.gz", optional: true))
    filtered_aln = tuple(meta, file("core-genome.filtered.aln.gz", optional: true))
    csv          = tuple(meta, file("panaroo/gene_presence_absence_roary.csv", optional: true))
    panaroo_csv  = tuple(meta, file("panaroo/gene_presence_absence.csv", optional: true))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
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
