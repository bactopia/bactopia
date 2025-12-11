/**
 * A tool to quickly download assemblies from NCBI's Assembly database.
 *
 * This process executes ncbigenomedownload to perform analysis
 *
 * @status stable
 * @keywords fasta, download, assembly
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-logic, resource-download
 * @citation ncbigenomedownload
 *
 * @input accessions
 * List of accessions (one per line) to download
 *
 * @output all        All downloaded files
 * @output gbk        GenBank format of the genomic sequence(s) in the assembly
 * @output fna        FASTA format of the genomic sequence(s) in the assembly.
 * @output rm         RepeatMasker output for eukaryotes.
 * @output features   Tab-delimited text file reporting locations and attributes for a subset of annotated features
 * @output gff        Annotation of the genomic sequence(s) in GFF3 format
 * @output faa        FASTA format of the accessioned protein products annotated on the genome assembly.
 * @output gpff       GenPept format of the accessioned protein products annotated on the genome assembly.
 * @output wgs_gbk    GenBank flat file format of the WGS master for the assembly
 * @output cds        FASTA format of the nucleotide sequences corresponding to all CDS features
 * @output rna        FASTA format of accessioned RNA products annotated on the genome assembly
 * @output rna_fna    FASTA format of the nucleotide sequences corresponding to all RNA features
 * @output report     Tab-delimited text file reporting the name, role and relationships of assembly units and sequences
 * @output stats      Tab-delimited text file reporting statistics for the assembly
 * @output accessions Accession list files
 * @output logs       Optional tool execution logs
 * @output nf_logs    Nextflow execution logs
 * @output versions   Software version information (YAML format)
 */
nextflow.preview.types = true

process NCBIGENOMEDOWNLOAD {
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    accessions : Path

    output:
    all        = files("*.gz")
    gbk        = tuple(meta, files("*_genomic.gbff.gz", optional: true))
    fna        = tuple(meta, files("*_genomic.fna.gz", optional: true))
    rm         = tuple(meta, files("*_rm.out.gz", optional: true))
    features   = tuple(meta, files("*_feature_table.txt.gz", optional: true))
    gff        = tuple(meta, files("*_genomic.gff.gz", optional: true))
    faa        = tuple(meta, files("*_protein.faa.gz", optional: true))
    gpff       = tuple(meta, files("*_protein.gpff.gz", optional: true))
    wgs_gbk    = tuple(meta, files("*_wgsmaster.gbff.gz", optional: true))
    cds        = tuple(meta, files("*_cds_from_genomic.fna.gz", optional: true))
    rna        = tuple(meta, files("*_rna.fna.gz", optional: true))
    rna_fna    = tuple(meta, files("*_rna_from_genomic.fna.gz", optional: true))
    report     = tuple(meta, files("*_assembly_report.txt", optional: true))
    stats      = tuple(meta, files("*_assembly_stats.txt", optional: true))
    accessions = tuple(meta, files("accession-*.txt", optional: true))
    logs       = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs   = tuple(meta, files(".command.*"))
    versions   = tuple(meta, file("versions.yml"))

    script:
    meta = [:]
    meta.id = task.ext.meta_id
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$(echo \$(ncbi-genome-download --version 2>&1) | sed 's/ncbi-genome-download //')
    END_VERSIONS
    """
}
