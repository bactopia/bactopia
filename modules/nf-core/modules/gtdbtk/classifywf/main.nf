// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

def VERSION = '1.5.0' // When using stubs for the GTDB database, the version info isn't printed.

process GTDBTK_CLASSIFYWF {
    tag "${meta.assembler}-${meta.id}"
    label 'process_high'
    publishDir "${publish_dir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:getSoftwareName(task.process, options.full_software_name), is_module: options.is_module, publish_to_base: options.publish_to_base) }

    conda (params.enable_conda ? "bioconda::gtdbtk=1.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gtdbtk:1.5.0--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/gtdbtk:1.5.0--pyhdfd78af_0"
    }

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

    stub:
    """
    touch gtdbtk.${meta.assembler}-${meta.id}.stub.summary.tsv
    touch gtdbtk.${meta.assembler}-${meta.id}.stub.classify.tree.gz
    touch gtdbtk.${meta.assembler}-${meta.id}.stub.markers_summary.tsv
    touch gtdbtk.${meta.assembler}-${meta.id}.stub.msa.fasta.gz
    touch gtdbtk.${meta.assembler}-${meta.id}.stub.user_msa.fasta
    touch gtdbtk.${meta.assembler}-${meta.id}.stub.filtered.tsv
    touch gtdbtk.${meta.assembler}-${meta.id}.log
    touch gtdbtk.${meta.assembler}-${meta.id}.warnings.log
    touch gtdbtk.${meta.assembler}-${meta.id}.failed_genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo "$VERSION")
    END_VERSIONS
    """
}
