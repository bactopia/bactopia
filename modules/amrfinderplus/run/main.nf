/**
 * Identify antimicrobial resistance and virulence genes in gene or protein sequences.
 *
 * Uses [AMRFinder+](https://github.com/ncbi/amr) to screen nucleotide or protein
 * sequences against NCBI's [Reference Gene Database](https://www.ncbi.nlm.nih.gov/pathogens/isolates#/refgene/).
 * It identifies AMR genes, resistance-associated point mutations, and select other classes of
 * genes using protein annotations and/or assembled nucleotide sequences.
 *
 * @status stable
 * @keywords bacteria, fasta, antimicrobial resistance, virulence, ncbi, amr, genes, proteins
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent
 * @citation amrfinderplus
 *
 * @note Requires external database to be available
 *
 * @input tuple(meta, genes, proteins, gff)
 * - `meta`: Groovy Map containing sample information
 * - `genes`: Nucleotide sequences of genes in FASTA format
 * - `proteins`: Amino acid sequences of proteins in FASTA format
 * - `gff`: Genome annotation in GFF3 format
 *
 * @input db
 * A compressed tarball of the AMRFinderPlus database to query
 *
 * @output report          Tab-delimited report of identified AMR genes and virulence factors
 * @output mutation_report Tab-delimited report of identified point mutations (if applicable)
 * @output logs            Optional software execution logs containing warnings/errors
 * @output nf_logs         Nextflow execution scripts and logs for debugging
 * @output versions        A YAML formatted file with software versions
 */
nextflow.preview.types = true

process AMRFINDERPLUS_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, genes, proteins, gff) : Tuple<Map, Set<Path>, Set<Path>, Set<Path>>
    db                            : Path

    output:
    report          = tuple(meta, files("${prefix}.tsv"))
    mutation_report = tuple(meta, files("${prefix}-mutations.tsv", optional: true))
    logs            = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs         = tuple(meta, files(".command.*"))
    versions        = tuple(meta, files("versions.yml"))

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

    // WF specific parameters
    def fna_is_compressed = genes.toList()[0].getName().endsWith(".gz") ? true : false
    def faa_is_compressed = proteins.toList()[0].getName().endsWith(".gz") ? true : false
    def gff_is_compressed = gff.toList()[0].getName().endsWith(".gz") ? true : false
    organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    fna_name = genes.toList()[0].getName().replace(".gz", "")
    faa_name = proteins.toList()[0].getName().replace(".gz", "")
    gff_name = gff.toList()[0].getName().replace(".gz", "")
    annotation_format = gff_name.endsWith(".gff") ? "prokka" : "bakta"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${fna_is_compressed}" == "true" ]; then
        gzip -c -d ${genes} > ${fna_name}
    fi

    if [ "${faa_is_compressed}" == "true" ]; then
        gzip -c -d ${proteins} > ${faa_name}
    fi

    if [ "${gff_is_compressed}" == "true" ]; then
        gzip -c -d ${gff} > ${gff_name}
    fi

    # Extract database
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        AMRFINDER_DB=\$(find database/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    else
        AMRFINDER_DB=\$(find ${db}/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    fi

    # Full AMRFinderPlus search combining results
    amrfinder \\
        --nucleotide ${fna_name} \\
        --protein ${faa_name} \\
        --gff ${gff_name} \\
        --annotation_format ${annotation_format} \\
        ${organism_param} \\
        ${task.ext.args} \\
        --database \$AMRFINDER_DB \\
        --threads ${task.cpus} \\
        --name ${prefix} > ${prefix}.tsv

    # Clean up
    DB_VERSION=\$(echo \$(echo \$(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    rm -rf database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$DB_VERSION
    END_VERSIONS
    """
}
