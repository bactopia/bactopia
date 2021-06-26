nextflow.enable.dsl = 2

process CALL_VARIANTS_AUTO {
    /*
    Identify variants (SNPs/InDels) against one or more reference genomes selected based
    on their Mash distance from the input.
    */
    tag "${sample} - ${reference_name}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/variants/auto", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${reference_name}/*"

    input:
    tuple val(sample), val(single_end), path(fq), path(reference)

    output:
    path "${reference_name}/*"
    path "${task.process}/*" optional true

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getSimpleName().split("${sample}-")[1].split(/\./)[0]
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    '''
    LOG_DIR="!{task.process}"
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
    echo "# Snippy Version" >> ${LOG_DIR}/!{task.process}.versions
    snippy --version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    # Print captured STDERR incase of exit
    function print_stderr {
        cat .command.err 1>&2
        ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
    }
    trap print_stderr EXIT

    # Verify AWS files were staged
    if [[ ! -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "true" ]; then
            check-staging.py --fq1 !{fq[0]} --extra !{reference} --is_single
        else
            check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{reference}
        fi
    fi

    snippy !{fastq} \
        --ref !{reference} \
        --cpus !{task.cpus} \
        --ram !{snippy_ram} \
        --outdir !{reference_name} \
        --prefix !{sample} \
        --mapqual !{params.mapqual} \
        --basequal !{params.basequal} \
        --mincov !{params.mincov} \
        --minfrac !{params.minfrac} \
        --minqual !{params.minqual} \
        --maxsoft !{params.maxsoft} !{bwaopt} !{fbopt} > ${LOG_DIR}/snippy.out 2> ${LOG_DIR}/snippy.err

    # Add GenBank annotations
    echo "# vcf-annotator Version" >> ${LOG_DIR}/!{task.process}.versions
    vcf-annotator --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    vcf-annotator !{reference_name}/!{sample}.vcf !{reference} > !{reference_name}/!{sample}.annotated.vcf 2> ${LOG_DIR}/vcf-annotator.err

    # Get per-base coverage
    echo "# bedtools Version" >> ${LOG_DIR}/!{task.process}.versions
    bedtools --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    grep "^##contig" !{reference_name}/!{sample}.vcf > !{reference_name}/!{sample}.full-coverage.txt
    genomeCoverageBed -ibam !{reference_name}/!{sample}.bam -d >> !{reference_name}/!{sample}.full-coverage.txt 2> ${LOG_DIR}/genomeCoverageBed.err
    cleanup-coverage.py !{reference_name}/!{sample}.full-coverage.txt > !{reference_name}/!{sample}.coverage.txt
    rm !{reference_name}/!{sample}.full-coverage.txt

    echo "here 6"
    # Mask low coverage regions
    mask-consensus.py !{sample} !{reference_name} \
                    !{reference_name}/!{sample}.consensus.subs.fa \
                    !{reference_name}/!{sample}.subs.vcf \
                    !{reference_name}/!{sample}.coverage.txt \
                    --mincov !{params.mincov}
    echo "here 7"
    # Clean Up
    rm -rf !{reference_name}/reference !{reference_name}/ref.fa* !{reference_name}/!{sample}.vcf.gz*
    echo "here 8"
    if [[ !{params.compress} == "true" ]]; then
        find !{reference_name}/ -type f -not -name "*.bam*" -and -not -name "*.log*" -and -not -name "*.txt*" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
        pigz -n --best -p !{task.cpus} !{reference_name}/!{sample}.coverage.txt
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
        cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    reference_name = "ref_name"
    """
    echo True
    mkdir ${reference_name}
    mkdir ${task.process}
    touch ${reference_name}/*
    touch ${task.process}/*
    """
}


//###############
//Module testing
//###############

workflow test {
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        path(params.fq),
        path(params.reference)
        ])
    call_variants_auto(TEST_PARAMS_CH)
}
