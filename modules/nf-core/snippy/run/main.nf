// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'snippy')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
publish_dir = options.is_module ? "${publish_dir}/snippy" : publish_dir
conda_tools = "bioconda::snippy=4.6.0 bioconda::vcf-annotator=0.7 bioconda::rename=1.601 bioconda::snpeff=5.0 conda-forge::openjdk=11.0.13"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process SNIPPY_RUN {
    tag "${meta.id} - ${reference_name}"
    label "process_low"
    label "base_mem_4gb"
    label "max_cpu_75"

    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container 'quay.io/bactopia/call_variants:2.2.0'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("results/${prefix}.aligned.fa.gz")                                    , emit: aligned_fa
    tuple val(meta), path("results/${prefix}.annotated.vcf.gz")                                 , emit: annotated_vcf
    tuple val(meta), path("results/${prefix}.bam")                                              , emit: bam
    tuple val(meta), path("results/${prefix}.bam.bai")                                          , emit: bai
    tuple val(meta), path("results/${prefix}.bed.gz")                                           , emit: bed
    tuple val(meta), path("results/${prefix}.consensus.fa.gz")                                  , emit: consensus_fa
    tuple val(meta), path("results/${prefix}.consensus.subs.fa.gz")                             , emit: consensus_subs_fa
    tuple val(meta), path("results/${prefix}.consensus.subs.masked.fa.gz")                      , emit: consensus_subs_masked_fa
    tuple val(meta), path("results/${prefix}.coverage.txt.gz")                                  , emit: coverage
    tuple val(meta), path("results/${prefix}.csv.gz")                                           , emit: csv
    tuple val(meta), path("results/${prefix}.filt.vcf.gz")                                      , emit: filt_vcf
    tuple val(meta), path("results/${prefix}.gff.gz")                                           , emit: gff
    tuple val(meta), path("results/${prefix}.html")                                             , emit: html
    tuple val(meta), path("results/${prefix}.raw.vcf.gz")                                       , emit: raw_vcf
    tuple val(meta), path("results/${prefix}.subs.vcf.gz")                                      , emit: subs_vcf
    tuple val(meta), path("results/${prefix}.tab")                                              , emit: tab
    tuple val(meta), path("results/${prefix}.txt")                                              , emit: txt
    tuple val(meta), path("results/${prefix}.vcf.gz")                                           , emit: vcf
    path "*.{log,err}"                                                          , optional: true, emit: logs
    path ".command.*"                                                                           , emit: nf_logs
    path "versions.yml"                                                                         , emit: versions

    when:
    meta.runtype != "ont"

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def read_inputs = meta.single_end ? "--se ${reads[0]}" : "--R1 ${reads[0]} --R2 ${reads[1]}"
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    def final_reference = reference.getName().replace(".gz", "")
    reference_name = reference.getSimpleName()
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $reference > $final_reference
    fi

    snippy \\
        $options.args \\
        --cpus $task.cpus \\
        --outdir $prefix \\
        --reference $final_reference \\
        --prefix $prefix \\
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
    rm -rf ${prefix}/reference ${prefix}/ref.fa* ${prefix}/${prefix}.vcf.gz*

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
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/snippy //')
        vcf-annotator: \$(echo \$(vcf-annotator --version 2>&1) | sed 's/vcf-annotator.py //')
    END_VERSIONS
    """
}
