/**
 * Determine the phylogroup of Escherichia coli isolates.
 *
 * Uses [ClermonTyping](https://github.com/A-BN/ClermonTyping) to perform in silico PCR
 * detection of specific marker genes (arpA, chuA, yjaA, TspE4.C2). This assigns the
 * isolate to one of the main *E. coli* phylogroups (A, B1, B2, C, D, E, F, G, or Cryptic).
 *
 * @status stable
 * @keywords bacteria, escherichia coli, typing, phylotyping, pcr, phylogroup
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation clermontyping
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv           A tab-delimited report containing the assigned phylogroup and genotype
 * @output supplemental  Directory containing detailed analysis files and HTML reports
 * @output logs          Optional software execution logs containing warnings/errors
 * @output nf_logs       Nextflow execution scripts and logs for debugging
 * @output versions      A YAML formatted file with software versions
 */
nextflow.preview.types = true

process CLERMONTYPING {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

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
    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    clermonTyping.sh \\
        --fasta ${assembly_name} \\
        --name supplemental \\
        ${task.ext.args}

    # Remove temporary files and rename outputs
    rm supplemental/${assembly_name} ${assembly_name}
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
