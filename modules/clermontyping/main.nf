/**
 * In silico PCR for typing Escherichia coli isolates.
 *
 * This process executes clermontyping to perform analysis
 *
 * @status stable
 * @keywords escherichia coli, phylotyping, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output, compression, conditional-logic
 * @citation clermontyping
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Genome assembly in FASTA format
 *
 * @output tsv          ClermonTyping phylogroup results
 * @output supplemental Supplemental
 * @output logs         Optional tool execution logs
 * @output nf_logs      Nextflow execution logs
 * @output versions     Software version information (YAML format)
 */
nextflow.preview.types = true

process CLERMONTYPING {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    output:
    tsv          = tuple(meta, file("${prefix}.tsv"))
    supplemental = tuple(meta, files("supplemental/*"))
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
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    clermonTyping.sh \\
        --fasta ${fasta_name} \\
        --name supplemental \\
        ${task.ext.args}

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
