/**
 * Download genome assemblies and annotation files from NCBI Datasets.
 *
 * Uses [genome-dl](https://github.com/rpetit3/genome-dl) to query the NCBI Datasets v2 REST API
 * for assembly metadata, then download the requested file formats directly from the NCBI FTP
 * site. Accessions are resolved to their latest version before download, and species queries
 * are subsampled to the first `--limit` assemblies in NCBI relevance order (reference first).
 *
 * @status stable
 * @keywords ncbi, datasets, download, genome, assembly, fasta, utility
 * @tags complexity:moderate input-type:single output-type:multiple features:internet-access,resource-download,conditional-logic
 * @citation genome_dl
 *
 * @note Internet Required
 * Queries the NCBI Datasets v2 REST API and downloads from the NCBI FTP site. Setting the
 * `NCBI_API_KEY` environment variable raises the API rate limit from 5 to 10 requests per second.
 *
 * @note Species Downloads Are Capped
 * `--limit` defaults to 100 so a broad `--species` cannot accidentally pull the tens of
 * thousands of assemblies NCBI holds for common taxa. Raise it, or use `--limit 0` for no
 * limit, only when that is genuinely intended.
 *
 * @input accessions?
 * A path to a text file containing a list of NCBI Assembly accessions (one per line)
 *
 * @output record(meta, fna?, gbff?, wgs_gbk?, gff?, gtf?, faa?, gpff?, cds?, translated_cds?, rna?, features?, report?, stats?, metadata?, summary?, json?, results, logs, nf_logs, versions)
 * - `fna?`: FASTA format of the genomic nucleotide sequence(s) (*.fna.gz)
 * - `gbff?`: GenBank format of the genomic sequence(s) (*.gbff.gz)
 * - `wgs_gbk?`: GenBank flat file format of the WGS master (*.wgsmaster.gbff.gz)
 * - `gff?`: Annotation of the genomic sequence(s) in GFF3 format (*.gff.gz)
 * - `gtf?`: Annotation of the genomic sequence(s) in GTF format (*.gtf.gz)
 * - `faa?`: FASTA format of the accessioned protein products (*.faa.gz)
 * - `gpff?`: GenPept format of the accessioned protein products (*.gpff.gz)
 * - `cds?`: FASTA format of the nucleotide sequences corresponding to all CDS features
 * - `translated_cds?`: FASTA format of the protein sequences corresponding to all CDS features
 * - `rna?`: FASTA format of the nucleotide sequences corresponding to all RNA features
 * - `features?`: Tab-delimited text file reporting locations and attributes for a subset of features
 * - `report?`: Tab-delimited text file reporting assembly unit names, roles, and relationships
 * - `stats?`: Tab-delimited text file reporting assembly statistics
 * - `metadata?`: Tab-delimited NCBI Datasets metadata for each downloaded assembly
 * - `summary?`: Human-readable run summary of the version, parameters, and results
 * - `json?`: Machine-readable run report of the parameters, results, and per-assembly metadata
 */
nextflow.enable.types = true

// bactopia-lint: ignore M017,M026
process GENOMEDL {
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    accessions : Path?

    stage:
    stageAs accessions, 'staging/accessions/*'

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        // Use the [0-9] to separate a genomic FASTA from the cds.fna.gz/rna.fna.gz variants
        fna: files("*[0-9].fna.gz", optional: true),
        gbff: files("*[0-9].gbff.gz", optional: true),
        wgs_gbk: files("*.wgsmaster.gbff.gz", optional: true),
        gff: files("*.gff.gz", optional: true),
        gtf: files("*.gtf.gz", optional: true),
        faa: files("*[0-9].faa.gz", optional: true),
        gpff: files("*.gpff.gz", optional: true),
        cds: files("*.cds.fna.gz", optional: true),
        translated_cds: files("*.translated_cds.faa.gz", optional: true),
        rna: files("*.rna.fna.gz", optional: true),
        features: files("*.feature_table.txt.gz", optional: true),
        report: files("*.assembly_report.txt", optional: true),
        stats: files("*.assembly_stats.txt", optional: true),
        metadata: files("*-metadata.tsv", optional: true),
        summary: files("*-summary.txt", optional: true),
        json: files("*.json", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("*.gz", optional: true),
            files("*.txt", optional: true),
            files("*.tsv", optional: true),
            files("*.json", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: task.ext.meta_id
    meta = record(
        id: task.ext.meta_id,
        name: task.ext.meta_id,
        limit: task.ext.meta_limit,
        accession: task.ext.meta_accession,
        species: task.ext.meta_species,
        scope: task.ext.scope,
        process_name: task.ext.process_name,
        output_dir: task.ext.process_name,
        logs_dir: "${task.ext.process_name}/logs"
    )

    def has_accession = task.ext.meta_accession != null
    def has_accessions = accessions != null
    def has_species = task.ext.meta_species != null
    def opts = "${task.ext.args} --outdir ./ --cpus ${task.cpus} --max-attempts ${task.ext.max_retry}"
    """
    # genome-dl can't mix --accession, --accessions, and --species, so run for each
    if [ "${has_accession}" == "true" ]; then
        genome-dl ${opts} \\
            --prefix ${prefix}-accession \\
            --accession ${task.ext.meta_accession}
    fi

    if [ "${has_accessions}" == "true" ]; then
        genome-dl ${opts} \\
            --prefix ${prefix}-accessions \\
            --accessions ${accessions}
    fi

    if [ "${has_species}" == "true" ]; then
        genome-dl ${opts} ${task.ext.args2} \\
            --prefix ${prefix}-species \\
            --species "${task.ext.meta_species}"
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomedl: \$(echo \$(genome-dl --version 2>&1) | sed 's/.*version //')
    END_VERSIONS
    """
}
