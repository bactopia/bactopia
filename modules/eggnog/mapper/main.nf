/**
 * Functional annotation of proteins using eggNOG orthology data.
 *
 * This process executes eggnog_mapper to perform analysis
 *
 * @status stable
 * @keywords eggnog, annotation, orthology, proteins
 * @tags complexity:complex input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, database-dependent
 * @citation eggnog_mapper
 *
 * @note Requires external database to be available
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Protein sequences in FASTA format
 *
 * @input db
 * eggNOG database directory or tarball
 *
 * @output hits           Diamond/MMseqs2 hits file
 * @output seed_orthologs Seed orthologs file
 * @output annotations    Tab-delimited annotations file
 * @output xlsx           Excel format annotations (optional)
 * @output orthologs      Orthologs file (optional)
 * @output genepred       Predicted genes (optional)
 * @output gff            GFF format annotations (optional)
 * @output no_anno        Sequences without annotations (optional)
 * @output pfam           PFAM annotations (optional)
 * @output logs           Optional tool execution logs
 * @output nf_logs        Nextflow execution logs
 * @output versions       Software version information (YAML format)
 */
nextflow.preview.types = true

process EGGNOG_MAPPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>
    db             : Path

    output:
    hits           = tuple(meta, files("*.emapper.hits"))
    seed_orthologs = tuple(meta, files("*.emapper.seed_orthologs"))
    annotations    = tuple(meta, files("*.emapper.annotations"))
    xlsx           = tuple(meta, files("*.emapper.annotations.xlsx", optional: true))
    orthologs      = tuple(meta, files("*.emapper.orthologs", optional: true))
    genepred       = tuple(meta, files("*.emapper.genepred.fasta", optional: true))
    gff            = tuple(meta, files("*.emapper.gff", optional: true))
    no_anno        = tuple(meta, files("*.emapper.no_annotations.fasta", optional: true))
    pfam           = tuple(meta, files("*.emapper.pfam", optional: true))
    logs           = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs        = tuple(meta, files(".command.*"))
    versions       = tuple(meta, file("versions.yml"))

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
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        EGGNOG_DB=\$(find database/ -name "eggnog.db" | sed 's=eggnog.db==')
    else
        EGGNOG_DB=\$(find ${db}/ -name "eggnog.db" | sed 's=eggnog.db==')
    fi

    emapper.py \\
        ${task.ext.args} \\
        --cpu ${task.cpus} \\
        --data_dir \$EGGNOG_DB \\
        --output ${prefix} \\
        -i ${fasta}

    # Cleanup
    if [ "${is_tarball}" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
