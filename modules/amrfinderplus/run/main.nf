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
 * - `proteins`: Optional amino acid sequences of proteins in FASTA format
 * - `gff`: Optional genome annotation in GFF3 format
 *
 * @input db
 * A compressed tarball of the AMRFinderPlus database to query
 *
 * @output record(meta, report, mutation_report, results, logs, nf_logs, versions)
 * - `report`: A tab-delimited report of identified AMR genes and virulence factors
 * - `mutation_report`: Organism-specific point mutations associated with antimicrobial resistance
 */
nextflow.preview.types = true

process AMRFINDERPLUS_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, genes: Path, proteins: Path, gff: Path): Record
    db: Path

    stage:
    stageAs 'genes/*', genes
    stageAs 'proteins/*', proteins
    stageAs 'gff/*', gff

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        report: file("${prefix}.tsv"),
        mutation_report: file("${prefix}-mutations.tsv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}-mutations.tsv", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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

    // Check for optional inputs
    def has_proteins = proteins != null
    def has_gff = gff != null

    // WF specific parameters
    def fna_cat = genes.getName().endsWith(".gz") ? "zcat" : "cat"
    def faa_cat = has_proteins ? proteins.getName().endsWith(".gz") ? "zcat" : "cat" : ""
    def gff_cat = has_gff ? gff.getName().endsWith(".gz") ? "zcat" : "cat" : ""
    organism_param = _meta.containsKey("organism") ? "--organism ${_meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    fna_name = "${prefix}.fna"
    faa_name = has_proteins ? "${prefix}.faa" : ""
    gff_name = has_gff ? "${prefix}.gff" : ""
    annotation_format = has_gff && gff_name.endsWith(".gff") ? "prokka" : "bakta"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false

    // Build optional parameters
    def protein_param = has_proteins ? "--protein ${faa_name}" : ""
    def gff_param = has_proteins && has_gff ? "--gff ${gff_name} --annotation_format ${annotation_format}" : ""
    """
    # Prepare input files
    ${fna_cat} ${genes} > ${fna_name}

    if [ "${has_proteins}" == "true" ]; then
        ${faa_cat} ${proteins} > ${faa_name}
    fi

    if [ "${has_gff}" == "true" ]; then
        ${gff_cat} ${gff} > ${gff_name}
    fi

    # Extract database
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        AMRFINDER_DB=\$(find database/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    else
        AMRFINDER_DB=\$(find ${db}/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    fi
    echo "Using AMRFINDER_DB: \$AMRFINDER_DB"

    # AMRFinderPlus search (with optional protein/gff inputs)
    amrfinder \\
        --nucleotide ${fna_name} \\
        ${protein_param} \\
        ${gff_param} \\
        ${organism_param} \\
        ${task.ext.args} \\
        --database \$AMRFINDER_DB \\
        --threads ${task.cpus} \\
        --name ${prefix} > ${prefix}.tsv

    # Clean up
    DB_VERSION=\$(echo \$(echo \$(amrfinder --database \$AMRFINDER_DB --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    rm -rf database/
    rm -rf ${fna_name} ${faa_name} ${gff_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$DB_VERSION
    END_VERSIONS
    """
}
