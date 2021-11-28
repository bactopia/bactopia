// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'

options     = initOptions(params.options ? params.options : [:], 'ncbigenomedownload')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process NCBIGENOMEDOWNLOAD {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::ncbi-genome-download=0.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.0--pyh864c0ab_1' :
        'quay.io/biocontainers/ncbi-genome-download:0.3.0--pyh864c0ab_1' }"

    input:
    val meta
    path accessions

    output:
    path("*.gz")                            , emit: all     , optional: true
    path("*_genomic.gbff.gz")               , emit: gbk     , optional: true
    path("*_genomic.fna.gz")                , emit: fna     , optional: true
    path("*_rm.out.gz")                     , emit: rm      , optional: true
    path("*_feature_table.txt.gz")          , emit: features, optional: true
    path("*_genomic.gff.gz")                , emit: gff     , optional: true
    path("*_protein.faa.gz")                , emit: faa     , optional: true
    path("*_protein.gpff.gz")               , emit: gpff    , optional: true
    path("*_wgsmaster.gbff.gz")             , emit: wgs_gbk , optional: true
    path("*_cds_from_genomic.fna.gz")       , emit: cds     , optional: true
    path("*_rna.fna.gz")                    , emit: rna     , optional: true
    path("*_rna_from_genomic.fna.gz")       , emit: rna_fna , optional: true
    path("*_assembly_report.txt")           , emit: report  , optional: true
    path("*_assembly_stats.txt")            , emit: stats   , optional: true
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs    , optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def has_accessions = accessions ? true : false
    def opts = "${options.args} --output-folder ./ --flat-output -p ${task.cpus} -r 50"
    """
    if [ "${meta.species}" != "null" ]; then
        if [ "${meta.limit}" != "null" ]; then
            ncbi-genome-download $opts -g "${meta.species}" --dry-run | grep -v "Considering" > accession-list.txt
            shuf accession-list.txt | head -n ${meta.limit} | cut -f 1,1  > accession-subset.txt
            ncbi-genome-download $opts -A accession-subset.txt
        else
            ncbi-genome-download $opts -g "${meta.species}"
        fi
    fi

    if [ "${meta.accession}" != "null" ]; then
        ncbi-genome-download $opts -A ${meta.accession}
    fi

    if [ "${meta.has_accessions}" == "true" ]; then
        ncbi-genome-download $opts -A ${accessions}
    fi

    cat <<-END_VERSIONS > versions.yml
    ncbigenomedownload:
        ncbi-genome-download: \$( ncbi-genome-download --version )
    END_VERSIONS
    """
}
