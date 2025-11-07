process SNIPPY_RUN {
    tag "${prefix} - ${reference_name}"
    label "process_low"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(reads)
    tuple val(ref_meta), path(reference)

    output:
    tuple val(meta), path("${prefix}.aligned.fa.gz")              , emit: aligned_fa, optional: true
    tuple val(meta), path("${prefix}.vcf.gz")                     , emit: vcf, optional: true
    tuple val(meta), path("${prefix}.error.aligned.fa.gz")        , emit: aligned_fa_error, optional: true
    tuple val(meta), path("${prefix}.error.vcf.gz")               , emit: vcf_error, optional: true
    tuple val(meta), path("${prefix}.error.txt")                  , emit: error, optional: true
    tuple val(meta), path("${prefix}.annotated.vcf.gz")           , emit: annotated_vcf
    tuple val(meta), path("${prefix}.bam")                        , emit: bam, optional: true
    tuple val(meta), path("${prefix}.bam.bai")                    , emit: bai, optional: true
    tuple val(meta), path("${prefix}.bed.gz")                     , emit: bed
    tuple val(meta), path("${prefix}.consensus.fa.gz")            , emit: consensus_fa
    tuple val(meta), path("${prefix}.consensus.subs.fa.gz")       , emit: consensus_subs_fa
    tuple val(meta), path("${prefix}.consensus.subs.masked.fa.gz"), emit: consensus_subs_masked_fa
    tuple val(meta), path("${prefix}.coverage.txt.gz")            , emit: coverage
    tuple val(meta), path("${prefix}.csv.gz")                     , emit: csv
    tuple val(meta), path("${prefix}.filt.vcf.gz")                , emit: filt_vcf
    tuple val(meta), path("${prefix}.gff.gz")                     , emit: gff
    tuple val(meta), path("${prefix}.html")                       , emit: html
    tuple val(meta), path("${prefix}.raw.vcf.gz")                 , emit: raw_vcf
    tuple val(meta), path("${prefix}.subs.vcf.gz")                , emit: subs_vcf
    tuple val(meta), path("${prefix}.tab")                        , emit: tab
    tuple val(meta), path("${prefix}.txt")                        , emit: txt
    tuple val(meta), path("*.{log,err}")                          , emit: logs, optional: true
    tuple val(meta), path(".command.begin")                       , emit: nf_begin
    tuple val(meta), path(".command.err")                         , emit: nf_err
    tuple val(meta), path(".command.log")                         , emit: nf_log
    tuple val(meta), path(".command.out")                         , emit: nf_out
    tuple val(meta), path(".command.run")                         , emit: nf_run
    tuple val(meta), path(".command.sh")                          , emit: nf_sh
    tuple val(meta), path(".command.trace")                       , emit: nf_trace
    tuple val(meta), path("versions.yml")                         , emit: versions

    script:
    reference_name = reference.getSimpleName()
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${reference_name}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${reference_name}/logs"
    meta.process_name = task.ext.process_name
    def read_inputs = _meta.single_end ? "--se ${reads[0]}" : "--R1 ${reads[0]} --R2 ${reads[1]}"
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    def final_reference = reference.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $reference > $final_reference
    fi

    if ! head -n 1 $final_reference | grep "^LOCUS"; then
        echo "ERROR: Reference file (${reference}) does not appear to be a GenBank file"
        exit 1
    fi

    mkdir tmp_snippy/
    snippy \\
        $task.ext.args \\
        --cpus $task.cpus \\
        --outdir $prefix \\
        --reference $final_reference \\
        --prefix $prefix \\
        --tmpdir tmp_snippy/ \\
        $read_inputs

    # Add GenBank annotations
    vcf-annotator ${prefix}/${prefix}.vcf ${final_reference} > ${prefix}/${prefix}.annotated.vcf

    # Get per-base coverage
    grep "^##contig" ${prefix}/${prefix}.vcf > ${prefix}/${prefix}.full-coverage.txt
    genomeCoverageBed -ibam ${prefix}/${prefix}.bam -d >> ${prefix}/${prefix}.full-coverage.txt
    cleanup-coverage.py ${prefix}/${prefix}.full-coverage.txt > ${prefix}/${prefix}.coverage.txt
    rm ${prefix}/${prefix}.full-coverage.txt

    # Mask low coverage regions
    mask-consensus.py \\
        ${prefix} \\
        ${reference_name} \\
        ${prefix}/${prefix}.consensus.subs.fa \\
        ${prefix}/${prefix}.subs.vcf \\
        ${prefix}/${prefix}.coverage.txt \\
        --mincov ${task.ext.mincov} > ${prefix}/${prefix}.consensus.subs.masked.fa

    # Clean Up
    rm -rf tmp_snippy/ ${prefix}/reference ${prefix}/ref.fa* ${prefix}/${prefix}.vcf.gz* ${final_reference}

    if [[ ${task.ext.snippy_remove_bam} == "true" ]]; then
        rm ${prefix}/${prefix}.bam ${prefix}/${prefix}.bam.bai
    fi

    if [[ ${task.ext.skip_compression} == "false" ]]; then
        find ${prefix}/ -type f | \
            grep -v -E "\\.bam\$|\\.bai\$|\\.log\$|\\.txt\$|\\.html\$|\\.tab\$" | \
            xargs -I {} pigz -n --best -p ${task.cpus} {}
        pigz -n --best -p ${task.cpus} ${prefix}/${prefix}.coverage.txt
    fi
    mv ${prefix}/* ./

    # Check for SNPs
    TOTAL_SNPS=\$(wc -l < ${prefix}.tab)
    echo "Total SNPS === \$TOTAL_SNPS"
    if [[ "\$TOTAL_SNPS" -eq 1 ]]; then
        # No SNPs were found
        mv ${prefix}.aligned.fa.gz ${prefix}.error.aligned.fa.gz
        mv ${prefix}.vcf.gz ${prefix}.error.vcf.gz
        echo "No SNPs found using reference ${final_reference}, downstream analysis is discontinued for ${prefix}" | \
        sed 's/^\\s*//' > ${prefix}-error.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/bedtools v//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy //')
        vcf-annotator: \$(echo \$(vcf-annotator --version 2>&1) | sed 's/vcf-annotator.py //')
    END_VERSIONS
    """
}
