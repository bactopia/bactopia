process SNIPPY_RUN {
    tag "${meta.id} - ${reference_name}"
    label "process_low"

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)
    tuple val(ref_meta), path(reference)

    output:
    tuple val(meta), path("results/${prefix}.aligned.fa.gz")                  , emit: aligned_fa
    tuple val(meta), path("results/${prefix}.annotated.vcf.gz")               , emit: annotated_vcf
    tuple val(meta), path("results/${prefix}.bam")                            , emit: bam, optional: true
    tuple val(meta), path("results/${prefix}.bam.bai")                        , emit: bai, optional: true
    tuple val(meta), path("results/${prefix}.bed.gz")                         , emit: bed
    tuple val(meta), path("results/${prefix}.consensus.fa.gz")                , emit: consensus_fa
    tuple val(meta), path("results/${prefix}.consensus.subs.fa.gz")           , emit: consensus_subs_fa
    tuple val(meta), path("results/${prefix}.consensus.subs.masked.fa.gz")    , emit: consensus_subs_masked_fa
    tuple val(meta), path("results/${prefix}.coverage.txt.gz")                , emit: coverage
    tuple val(meta), path("results/${prefix}.csv.gz")                         , emit: csv
    tuple val(meta), path("results/${prefix}.filt.vcf.gz")                    , emit: filt_vcf
    tuple val(meta), path("results/${prefix}.gff.gz")                         , emit: gff
    tuple val(meta), path("results/${prefix}.html")                           , emit: html
    tuple val(meta), path("results/${prefix}.raw.vcf.gz")                     , emit: raw_vcf
    tuple val(meta), path("results/${prefix}.subs.vcf.gz")                    , emit: subs_vcf
    tuple val(meta), path("results/${prefix}.tab")                            , emit: tab
    tuple val(meta), path("results/${prefix}.txt")                            , emit: txt
    tuple val(meta), path("results/${prefix}.vcf.gz")                         , emit: vcf
    tuple val(meta), path("*.{log,err}")                                      , emit: logs, optional: true
    tuple val(meta), path(".command.begin")                                   , emit: nf_begin
    tuple val(meta), path(".command.err")                                     , emit: nf_err
    tuple val(meta), path(".command.log")                                     , emit: nf_log
    tuple val(meta), path(".command.out")                                     , emit: nf_out
    tuple val(meta), path(".command.run")                                     , emit: nf_run
    tuple val(meta), path(".command.sh")                                      , emit: nf_sh
    tuple val(meta), path(".command.trace")                                   , emit: nf_trace
    tuple val(meta), path("versions.yml")                                     , emit: versions

    when:
    meta.runtype != "ont"

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    reference_name = reference.getSimpleName()
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${reference_name}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${reference_name}/logs"
    meta.process_name = task.ext.process_name
    def read_inputs = meta.single_end ? "--se ${reads[0]}" : "--R1 ${reads[0]} --R2 ${reads[1]}"
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
        --mincov ${params.mincov} > ${prefix}/${prefix}.consensus.subs.masked.fa

    # Clean Up
    rm -rf tmp_snippy/ ${prefix}/reference ${prefix}/ref.fa* ${prefix}/${prefix}.vcf.gz* ${final_reference}

    if [[ ${params.snippy_remove_bam} == "true" ]]; then
        rm ${prefix}/${prefix}.bam ${prefix}/${prefix}.bam.bai
    fi

    if [[ ${params.skip_compression} == "false" ]]; then
        find ${prefix}/ -type f | \
            grep -v -E "\\.bam\$|\\.bai\$|\\.log\$|\\.txt\$|\\.html\$|\\.tab\$" | \
            xargs -I {} pigz -n --best -p ${task.cpus} {}
        pigz -n --best -p ${task.cpus} ${prefix}/${prefix}.coverage.txt
    fi
    mv ${prefix}/ results/
    mv results/*.log ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/bedtools v//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy //')
        vcf-annotator: \$(echo \$(vcf-annotator --version 2>&1) | sed 's/vcf-annotator.py //')
    END_VERSIONS
    """
}
