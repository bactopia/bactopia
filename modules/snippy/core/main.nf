nextflow.preview.types = true

process SNIPPY_CORE {
    tag "${prefix}"
    label "process_medium"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, _vcf, _aligned_fa) : Tuple<Map, Path, Path>
    (_ref_meta, reference)   : Tuple<Map, Path>
    mask                     : Path

    output:
    supplemental   = tuple(meta, file("snippy-core/*"))
    aln            = tuple(meta, file("snippy-core/${prefix}.aln.gz"))
    full_aln       = tuple(meta, file("${prefix}.full.aln.gz"))
    clean_full_aln = tuple(meta, file("${prefix}-clean.full.aln.gz"))
    tab            = tuple(meta, file("snippy-core/${prefix}.tab.gz"))
    vcf            = tuple(meta, file("snippy-core/${prefix}.vcf.gz"))
    txt            = tuple(meta, file("snippy-core/${prefix}.txt"))
    samples        = tuple(meta, file("${reference_name}.samples.txt"))
    logs           = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin       = tuple(meta, file(".command.begin"))
    nf_err         = tuple(meta, file(".command.err"))
    nf_log         = tuple(meta, file(".command.log"))
    nf_out         = tuple(meta, file(".command.out"))
    nf_run         = tuple(meta, file(".command.run"))
    nf_sh          = tuple(meta, file(".command.sh"))
    nf_trace       = tuple(meta, file(".command.trace"))
    versions       = tuple(meta, file("versions.yml"))

    script:
    reference_name = reference.getSimpleName()
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = ""
    meta.logs_dir = "${meta.process_name}/logs"
    def mask_opt = mask ? "--mask ${mask[0]}" : ""
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    def final_reference = reference.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${reference} > ${final_reference}
    fi

    # Collect samples into necessary folders
    mkdir samples
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//;s=./==' > ${reference_name}.samples.txt
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//' | xargs -I {} bash -c 'mkdir samples/{}'
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//' | xargs -I {} bash -c 'gzip -cdf {}.vcf.gz > samples/{}/{}.vcf'
    find . -name "*.aligned.fa.gz" | sed 's/.aligned.fa.gz//' | xargs -I {} bash -c 'gzip -cdf {}.aligned.fa.gz > samples/{}/{}.aligned.fa'

    # Run snippy-core
    snippy-core \\
        ${task.ext.args} \\
        --ref ${final_reference} \\
        --prefix ${prefix} \\
        ${mask_opt} \\
        samples/*

    # Cleanup the alignment
    snippy-clean_full_aln ${prefix}.full.aln > ${prefix}-clean.full.aln
    rm -rf *ref.fa ${final_reference} samples/

    # Compress outputs
    if [[ ${task.ext.skip_compression} == "false" ]]; then
        pigz -n --best -p ${task.cpus} ${prefix}.aln
        pigz -n --best -p ${task.cpus} ${prefix}.full.aln
        pigz -n --best -p ${task.cpus} ${prefix}-clean.full.aln
        pigz -n --best -p ${task.cpus} ${prefix}.tab
        pigz -n --best -p ${task.cpus} ${prefix}.vcf
    fi

    # Move outputs
    mkdir snippy-core
    mv ${prefix}* snippy-core/
    mv snippy-core/${prefix}-clean.full.aln.gz snippy-core/${prefix}.full.aln.gz ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy //')
    END_VERSIONS
    """
}
