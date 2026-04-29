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
 * @input record(meta, fna, faa, gff)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Nucleotide sequences of genes in FASTA format
 * - `faa`: Optional amino acid sequences of proteins in FASTA format
 * - `gff`: Optional genome annotation in GFF3 format
 *
 * @input db
 * A compressed tarball of the AMRFinderPlus database to query
 *
 * @output record(meta, report, mutation_report?, results, logs, nf_logs, versions)
 * - `report`: A tab-delimited report of identified AMR genes and virulence factors
 * - `mutation_report?`: Organism-specific point mutations associated with antimicrobial resistance
 */
nextflow.enable.types = true

process AMRFINDERPLUS_RUN {
    tag prefix
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path,
        faa: Path,
        gff: Path
    )
    db: Path

    stage:
    stageAs fna, "staging/fna/*"
    stageAs faa, "staging/faa/*"
    stageAs gff, "staging/gff/*"

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
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name
    )

    // WF specific parameters
    def fna_cat = fna.getName().endsWith(".gz") ? "zcat" : "cat"
    def faa_cat = faa.getName().endsWith(".gz") ? "zcat" : "cat"
    def gff_cat = gff.getName().endsWith(".gz") ? "zcat" : "cat"
    def organism_param = _meta.organism != null ? "--organism ${_meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    def annotation_format = gff.getName().endsWith(".gff.gz") || gff.getName().endsWith(".gff") ? "prokka" : "bakta"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    # Prepare input files
    ${fna_cat} ${fna} > ${prefix}.fna
    ${faa_cat} ${faa} > ${prefix}.faa
    ${gff_cat} ${gff} > ${prefix}.gff

    # Extract database
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        AMRFINDER_DB=\$(find database/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    else
        AMRFINDER_DB=\$(find ${db}/ -name "AMR.LIB" | sed 's=AMR.LIB==')
    fi
    echo "Using AMRFINDER_DB: \$AMRFINDER_DB"
    DB_VERSION=\$(echo \$(echo \$(amrfinder --database \$AMRFINDER_DB --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))

    # AMRFinderPlus search (with optional protein/gff inputs)
    amrfinder \\
        --nucleotide ${prefix}.fna \\
        --protein ${prefix}.faa \\
        --gff ${prefix}.gff \\
        --annotation_format "${annotation_format}" \\
        ${organism_param} \\
        ${task.ext.args} \\
        --database \$AMRFINDER_DB \\
        --threads ${task.cpus} \\
        --name ${prefix} > ${prefix}.tsv

    # Cleanup
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database/
    fi
    rm -rf ${prefix}.fna ${prefix}.faa ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$DB_VERSION
    END_VERSIONS
    """
}
