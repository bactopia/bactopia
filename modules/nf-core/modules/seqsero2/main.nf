
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process SEQSERO2 {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:getSoftwareName(task.process, options.full_software_name), is_module: options.is_module, publish_to_base: options.publish_to_base) }

    conda (params.enable_conda ? "bioconda::seqsero2=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqsero2:1.2.1--py_0"
    } else {
        container "quay.io/biocontainers/seqsero2:1.2.1--py_0"
    }

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("results/*_log.txt")   , emit: log
    tuple val(meta), path("results/*_result.tsv"), emit: tsv
    tuple val(meta), path("results/*_result.txt"), emit: txt
    path "*.{stdout.txt,stderr.txt,log,err}"     , emit: logs, optional: true
    path ".command.*"                            , emit: nf_logs
    path "versions.yml"                          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def is_compressed_fna = seqs[0].getName().endsWith("fna.gz") ? true : false
    def seq_name = is_compressed_fna ? seqs[0].getName().replace(".gz", "") : "${seqs}"
    """
    if [ "$is_compressed_fna" == "true" ]; then
        gzip -c -d ${seqs[0]} > $seq_name
    fi
    SeqSero2_package.py \\
        $options.args \\
        -d results/ \\
        -n $prefix \\
        -p $task.cpus \\
        -i $seq_name

    cat <<-END_VERSIONS > versions.yml
    seqsero2:
        seqsero2: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
    """
}
