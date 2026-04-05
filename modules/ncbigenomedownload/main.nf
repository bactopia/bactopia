/**
 * Download assemblies and annotation files from NCBI's Assembly database.
 *
 * Uses [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) to efficiently fetch
 * one or more complete genome assemblies and their associated annotation and report files from
 * the NCBI FTP site based on accession numbers, species name, or assembly ID.
 *
 * @status stable
 * @keywords ncbi, download, genome, assembly, fasta, genbank, utility
 * @tags complexity:moderate input-type:single output-type:multiple features:internet-access,resource-download,conditional-logic
 * @citation ncbigenomedownload
 *
 * @input accessions?
 * A path to a text file containing a list of NCBI Assembly accession numbers (one per line)
 *
 * @output record(meta, gbff?, fna?, rm?, features?, gff?, faa?, gpff?, wgs_gbk?, cds?, rna?, rna_fna?, report?, stats?, accessions?, results, logs, nf_logs, versions)
 * - `gbff?`: GenBank format of the genomic sequence(s) (*_genomic.gbff.gz)
 * - `fna?`: FASTA format of the genomic nucleotide sequence(s) (*_genomic.fna.gz)
 * - `rm?`: RepeatMasker output for eukaryotes
 * - `features?`: Tab-delimited text file reporting locations and attributes for a subset of features
 * - `gff?`: Annotation of the genomic sequence(s) in GFF3 format (*_genomic.gff.gz)
 * - `faa?`: FASTA format of the accessioned protein products (*_protein.faa.gz)
 * - `gpff?`: GenPept format of the accessioned protein products
 * - `wgs_gbk?`: GenBank flat file format of the WGS master
 * - `cds?`: FASTA format of the nucleotide sequences corresponding to all CDS features
 * - `rna?`: FASTA format of accessioned RNA products
 * - `rna_fna?`: FASTA format of the nucleotide sequences corresponding to all RNA features
 * - `report?`: Tab-delimited text file reporting assembly unit names, roles, and relationships
 * - `stats?`: Tab-delimited text file reporting assembly statistics
 * - `accessions?`: The generated accession list files
 */
nextflow.preview.types = true

// bactopia-lint: ignore M017,M026
process NCBIGENOMEDOWNLOAD {
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    accessions : Path?

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        gbff: files("*_genomic.gbff.gz", optional: true),
        fna: files("*_genomic.fna.gz", optional: true),
        rm: files("*_rm.out.gz", optional: true),
        features: files("*_feature_table.txt.gz", optional: true),
        gff: files("*_genomic.gff.gz", optional: true),
        faa: files("*_protein.faa.gz", optional: true),
        gpff: files("*_protein.gpff.gz", optional: true),
        wgs_gbk: files("*_wgsmaster.gbff.gz", optional: true),
        cds: files("*_cds_from_genomic.fna.gz", optional: true),
        rna: files("*_rna.fna.gz", optional: true),
        rna_fna: files("*_rna_from_genomic.fna.gz", optional: true),
        report: files("*_assembly_report.txt", optional: true),
        stats: files("*_assembly_stats.txt", optional: true),
        accessions: files("accession-*.txt", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("*.gz", optional: true),
            files("*.txt", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    meta = [:]
    meta.id = task.ext.meta_id
    meta.name = task.ext.meta_id
    meta.limit = task.ext.meta_limit
    meta.accession = task.ext.meta_accession
    meta.species = task.ext.meta_species
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = "${meta.process_name}"
    meta.logs_dir = "${meta.process_name}/logs"

    def has_accessions = accessions ? true : false
    def opts = "${task.ext.args} --output-folder ./ --flat-output -p ${task.cpus} -r ${task.ext.max_retry}"
    """
    if [ "${meta.species}" != "null" ]; then
        if [ "${meta.limit}" != "null" ]; then
            ncbi-genome-download ${opts} -g "${meta.species}" --dry-run | grep -v "Considering" > accession-list.txt
            shuf accession-list.txt | head -n ${meta.limit} | cut -f 1,1  > accession-subset.txt
            ncbi-genome-download ${opts} -u "https://ftp.ncbi.nlm.nih.gov/genomes" -A accession-subset.txt
        else
            ncbi-genome-download ${opts} -g "${meta.species}"
        fi
    fi

    if [ "${meta.accession}" != "null" ]; then
        ncbi-genome-download ${opts} -u "https://ftp.ncbi.nlm.nih.gov/genomes" -A ${meta.accession}
    fi

    if [ "${has_accessions}" == "true" ]; then
        ncbi-genome-download ${opts} -u "https://ftp.ncbi.nlm.nih.gov/genomes" -A ${accessions}
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$(echo \$(ncbi-genome-download --version 2>&1) | sed 's/ncbi-genome-download //')
    END_VERSIONS
    """
}
