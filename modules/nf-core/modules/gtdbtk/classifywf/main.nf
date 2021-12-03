// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'gtdb')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

def VERSION = '1.5.0' // When using stubs for the GTDB database, the version info isn't printed.

process GTDBTK_CLASSIFYWF {
    tag "${meta.assembler}-${meta.id}"
    label 'process_high'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::gtdbtk=1.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:1.5.0--pyhdfd78af_0' :
        'quay.io/biocontainers/gtdbtk:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path("bins/*")
    tuple val(db_name), path("database/*")

    output:
    path "gtdbtk.${meta.assembler}-${meta.id}.*.summary.tsv"        ,emit: summary
    path "gtdbtk.${meta.assembler}-${meta.id}.*.classify.tree.gz"   ,emit: tree
    path "gtdbtk.${meta.assembler}-${meta.id}.*.markers_summary.tsv",emit: markers
    path "gtdbtk.${meta.assembler}-${meta.id}.*.msa.fasta.gz"       ,emit: msa
    path "gtdbtk.${meta.assembler}-${meta.id}.*.user_msa.fasta"     ,emit: user_msa
    path "gtdbtk.${meta.assembler}-${meta.id}.*.filtered.tsv"       ,emit: filtered
    path "gtdbtk.${meta.assembler}-${meta.id}.failed_genomes.tsv"   ,emit: failed
    path "*.{stdout.txt,stderr.txt,log,err}"                        ,emit: logs, optional: true
    path ".command.*"                                               ,emit: nf_logs
    path "versions.yml"                                             ,emit: versions

    script:
    def pplacer_scratch = params.gtdbtk_pplacer_scratch ? "--scratch_dir pplacer_tmp" : ""
    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    if [ ${pplacer_scratch} != "" ] ; then
        mkdir pplacer_tmp
    fi

    gtdbtk classify_wf \\
        $options.args \\
        --genome_dir bins \\
        --prefix "gtdbtk.${meta.assembler}-${meta.id}" \\
        --out_dir "\${PWD}" \\
        --cpus $task.cpus \\
        --pplacer_cpus $params.gtdbtk_pplacer_cpus \\
        $pplacer_scratch \\
        --min_perc_aa $params.gtdbtk_min_perc_aa \\
        --min_af $params.gtdbtk_min_af

    gzip "gtdbtk.${meta.assembler}-${meta.id}".*.classify.tree "gtdbtk.${meta.assembler}-${meta.id}".*.msa.fasta
    mv gtdbtk.log "gtdbtk.${meta.assembler}-${meta.id}.log"
    mv gtdbtk.warnings.log "gtdbtk.${meta.assembler}-${meta.id}.warnings.log"

    cat <<-END_VERSIONS > versions.yml
    gtdbtk_classifywf:
        gtdb-tk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
