process NCBIGENOMEDOWNLOAD {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    val meta
    path accessions

    output:
    path("*.gz")                            , emit: all
    path("*_genomic.gbff.gz")               , emit: gbk       , optional: true
    path("*_genomic.fna.gz")                , emit: fna       , optional: true
    path("*_rm.out.gz")                     , emit: rm        , optional: true
    path("*_feature_table.txt.gz")          , emit: features  , optional: true
    path("*_genomic.gff.gz")                , emit: gff       , optional: true
    path("*_protein.faa.gz")                , emit: faa       , optional: true
    path("*_protein.gpff.gz")               , emit: gpff      , optional: true
    path("*_wgsmaster.gbff.gz")             , emit: wgs_gbk   , optional: true
    path("*_cds_from_genomic.fna.gz")       , emit: cds       , optional: true
    path("*_rna.fna.gz")                    , emit: rna       , optional: true
    path("*_rna_from_genomic.fna.gz")       , emit: rna_fna   , optional: true
    path("*_assembly_report.txt")           , emit: report    , optional: true
    path("*_assembly_stats.txt")            , emit: stats     , optional: true
    path "accession-*.txt"                  , emit: accessions, optional: true
    path "versions.yml"                     , emit: versions
    path ".command.begin"                   , emit: begin
    path ".command.err"                     , emit: err
    path ".command.log"                     , emit: log
    path ".command.out"                     , emit: out
    path ".command.run"                     , emit: run
    path ".command.sh"                      , emit: sh
    path ".command.trace"                   , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def has_accessions = accessions ? true : false
    def opts = "${args} --output-folder ./ --flat-output -p ${task.cpus} -r ${params.max_retry}"
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
