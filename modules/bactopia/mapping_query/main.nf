nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "mapping_query"

process MAPPING_QUERY {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */
    tag "${sample}"
    label "max_cpus"
    label "mapping_query"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "mapping/*"

    input:
    tuple val(sample), val(single_end), path(fq)
    path(query)

    output:
    file "mapping/*"
    file "${PROCESS_NAME}/*" optional true

    shell:
    bwa_mem_opts = params.bwa_mem_opts ? params.bwa_mem_opts : ""
    bwa_aln_opts = params.bwa_aln_opts ? params.bwa_aln_opts : ""
    bwa_samse_opts = params.bwa_samse_opts ? params.bwa_samse_opts : ""
    bwa_sampe_opts = params.bwa_sampe_opts ? params.bwa_sampe_opts : ""
    '''
    LOG_DIR="!{PROCESS_NAME}"
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions

    # Print captured STDERR incase of exit
    function print_stderr {
        cat .command.err 1>&2
        ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
    }
    trap print_stderr EXIT

    # Verify AWS files were staged
    if [[ ! -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "true" ]; then
            check-staging.py --fq1 !{fq[0]} --extra !{query} --is_single
        else
            check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{query}
        fi
    fi

    avg_len=`seqtk fqchk !{fq[0]} | head -n 1 | sed -r 's/.*avg_len: ([0-9]+).*;.*/\\1/'`
    ls !{query}/* | xargs -I {} grep -H "^>" {} | awk '{print $1}' | sed 's/:>/\\t/; s=.*/==; s/\\..*\\t/\\t/' > mapping.txt
    cat !{query}/* > multifasta.fa

    echo "# bwa Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    bwa 2>&1 | grep "Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
    
    bwa index multifasta.fa > ${LOG_DIR}/bwa-index.out 2> ${LOG_DIR}/bwa-index.err
    if [ "${avg_len}" -gt "70" ]; then
        bwa mem -M -t !{task.cpus} !{bwa_mem_opts} multifasta.fa !{fq} > bwa.sam
    else
        if [ "!{single_end}" == "true" ]; then
            bwa aln -f bwa.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[0]} > ${LOG_DIR}/bwa-aln.out 2> ${LOG_DIR}/bwa-aln.err
            bwa samse -n !{params.bwa_n} !{bwa_samse_opts} multifasta.fa  bwa.sai !{fq[0]} > bwa.sam 2> ${LOG_DIR}/bwa-samse.err
        else
            bwa aln -f r1.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[0]} > ${LOG_DIR}/bwa-aln.out 2> ${LOG_DIR}/bwa-aln.err
            bwa aln -f r2.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[1]} >> ${LOG_DIR}/bwa-aln.out 2>> ${LOG_DIR}/bwa-aln.err
            bwa sampe -n !{params.bwa_n} !{bwa_sampe_opts} multifasta.fa  r1.sai r2.sai !{fq[0]} !{fq[1]} > bwa.sam 2> ${LOG_DIR}/bwa-sampe.err
        fi
    fi

    # Write per-base coverage
    echo "# samtools Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    samtools 2>&1 | grep "Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
    samtools view -bS bwa.sam | samtools sort -o cov.bam - > ${LOG_DIR}/samtools.out 2> ${LOG_DIR}/samtools.err

    echo "# bedtools Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    bedtools --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
    genomeCoverageBed -ibam cov.bam -d > cov.txt 2> ${LOG_DIR}/genomeCoverageBed.err
    split-coverages.py mapping.txt cov.txt --outdir mapping

    if [[ !{params.compress} == "true" ]]; then
        pigz --best -n -p !{task.cpus} mapping/*.txt
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp mapping.txt ${LOG_DIR}/
        cp .command.err ${LOG_DIR}/!{PROCESS_NAME}.err
        cp .command.out ${LOG_DIR}/!{PROCESS_NAME}.out
        cp .command.sh ${LOG_DIR}/!{PROCESS_NAME}.sh || :
        cp .command.trace ${LOG_DIR}/!{PROCESS_NAME}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir ${PROCESS_NAME}
    mkdir mapping
    touch ${PROCESS_NAME}/${sample}
    touch mapping/${sample}
    """
}
