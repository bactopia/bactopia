process SNIPPY_CORE {
    tag "${prefix}"
    label "process_medium"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(vcf), path(aligned_fa)
    tuple val(ref_meta), path(reference)
    path mask

    output:
    tuple val(meta), path("results/*")                          , emit: results
    tuple val(meta), path("results/${prefix}.aln.gz")           , emit: aln
    tuple val(meta), path("results/${prefix}.full.aln.gz")      , emit: full_aln
    tuple val(meta), path("results/${prefix}-clean.full.aln.gz"), emit: clean_full_aln
    tuple val(meta), path("results/${prefix}.tab.gz")           , emit: tab
    tuple val(meta), path("results/${prefix}.vcf.gz")           , emit: vcf
    tuple val(meta), path("results/${prefix}.txt")              , emit: txt
    tuple val(meta), path("${reference_name}.samples.txt")     , emit: samples
    tuple val(meta), path("*.{log,err}")                        , emit: logs, optional: true
    tuple val(meta), path(".command.begin")                     , emit: nf_begin
    tuple val(meta), path(".command.err")                       , emit: nf_err
    tuple val(meta), path(".command.log")                       , emit: nf_log
    tuple val(meta), path(".command.out")                       , emit: nf_out
    tuple val(meta), path(".command.run")                       , emit: nf_run
    tuple val(meta), path(".command.sh")                        , emit: nf_sh
    tuple val(meta), path(".command.trace")                     , emit: nf_trace
    tuple val(meta), path("versions.yml")                       , emit: versions

    script:
    reference_name = reference.getSimpleName()
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${reference_name}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${reference_name}/logs"
    meta.process_name = task.ext.process_name
    def mask_opt = mask ? "--mask ${mask[0]}" : ""
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    def final_reference = reference.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $reference > $final_reference
    fi

    # Collect samples into necessary folders
    mkdir samples
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//;s=./==' > ${reference_name}.samples.txt
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//' | xargs -I {} bash -c 'mkdir samples/{}'
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//' | xargs -I {} bash -c 'gzip -cdf {}.vcf.gz > samples/{}/{}.vcf'
    find . -name "*.aligned.fa.gz" | sed 's/.aligned.fa.gz//' | xargs -I {} bash -c 'gzip -cdf {}.aligned.fa.gz > samples/{}/{}.aligned.fa'

    # Run snippy-core
    snippy-core \\
        $task.ext.args \\
        --ref $final_reference \\
        --prefix $prefix \\
        $mask_opt \\
        samples/*

    # Cleanup the alignment
    snippy-clean_full_aln ${prefix}.full.aln > ${prefix}-clean.full.aln
    rm -rf *ref.fa ${final_reference} samples/

    # Compress outputs
    if [[ ${params.skip_compression} == "false" ]]; then
        pigz -n --best -p ${task.cpus} ${prefix}.aln
        pigz -n --best -p ${task.cpus} ${prefix}.full.aln
        pigz -n --best -p ${task.cpus} ${prefix}-clean.full.aln
        pigz -n --best -p ${task.cpus} ${prefix}.tab
        pigz -n --best -p ${task.cpus} ${prefix}.vcf
    fi

    # Move outputs
    mkdir results
    mv ${prefix}* results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy //')
    END_VERSIONS
    """
}
