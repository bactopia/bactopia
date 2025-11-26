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
