process NCBIGENOMEDOWNLOAD {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    path accessions

    output:
    path "*.gz"                                       , emit: all
    tuple val(meta), path("*_genomic.gbff.gz")        , emit: gbk       , optional: true
    tuple val(meta), path("*_genomic.fna.gz")         , emit: fna       , optional: true
    tuple val(meta), path("*_rm.out.gz")              , emit: rm        , optional: true
    tuple val(meta), path("*_feature_table.txt.gz")   , emit: features  , optional: true
    tuple val(meta), path("*_genomic.gff.gz")         , emit: gff       , optional: true
    tuple val(meta), path("*_protein.faa.gz")         , emit: faa       , optional: true
    tuple val(meta), path("*_protein.gpff.gz")        , emit: gpff      , optional: true
    tuple val(meta), path("*_wgsmaster.gbff.gz")      , emit: wgs_gbk   , optional: true
    tuple val(meta), path("*_cds_from_genomic.fna.gz"), emit: cds       , optional: true
    tuple val(meta), path("*_rna.fna.gz")             , emit: rna       , optional: true
    tuple val(meta), path("*_rna_from_genomic.fna.gz"), emit: rna_fna   , optional: true
    tuple val(meta), path("*_assembly_report.txt")    , emit: report    , optional: true
    tuple val(meta), path("*_assembly_stats.txt")     , emit: stats     , optional: true
    tuple val(meta), path("accession-*.txt")          , emit: accessions, optional: true
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    meta = task.ext.meta
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = "${meta.process_name}"
    meta.logs_dir = "${meta.process_name}/logs"
    def has_accessions = accessions ? true : false
    def opts = "${task.ext.args} --output-folder ./ --flat-output -p ${task.cpus} -r ${task.ext.max_retry}"
    """
    if [ "${meta.species}" != "null" ]; then
        if [ "${meta.limit}" != "null" ]; then
            ncbi-genome-download $opts -g "${meta.species}" --dry-run | grep -v "Considering" > accession-list.txt
            shuf accession-list.txt | head -n ${meta.limit} | cut -f 1,1  > accession-subset.txt
            ncbi-genome-download $opts -u "https://ftp.ncbi.nlm.nih.gov/genomes" -A accession-subset.txt
        else
            ncbi-genome-download $opts -g "${meta.species}"
        fi
    fi

    if [ "${meta.accession}" != "null" ]; then
        ncbi-genome-download $opts -u "https://ftp.ncbi.nlm.nih.gov/genomes" -A ${meta.accession}
    fi

    if [ "$has_accessions" == "true" ]; then
        ncbi-genome-download $opts -u "https://ftp.ncbi.nlm.nih.gov/genomes" -A $accessions
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$(echo \$(ncbi-genome-download --version 2>&1) | sed 's/ncbi-genome-download //')
    END_VERSIONS
    """
}
