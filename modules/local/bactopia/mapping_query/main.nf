nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "mapping_query"

process MAPPING_QUERY {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */
    tag "${sample}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME) }

    input:
    tuple val(sample), val(single_end), path(fq)
    path(query)

    output:
    path "results/*"
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    bwa_mem_opts = params.bwa_mem_opts ? params.bwa_mem_opts : ""
    bwa_aln_opts = params.bwa_aln_opts ? params.bwa_aln_opts : ""
    bwa_samse_opts = params.bwa_samse_opts ? params.bwa_samse_opts : ""
    bwa_sampe_opts = params.bwa_sampe_opts ? params.bwa_sampe_opts : ""
    '''
    avg_len=`seqtk fqchk !{fq[0]} | head -n 1 | sed -r 's/.*avg_len: ([0-9]+).*;.*/\\1/'`
    ls !{query}/* | xargs -I {} grep -H "^>" {} | awk '{print $1}' | sed 's/:>/\\t/; s=.*/==; s/\\..*\\t/\\t/' > mapping.txt
    cat !{query}/* > multifasta.fa

    bwa index multifasta.fa > bwa-index.stdout.txt 2> bwa-index.stderr.txt
    if [ "${avg_len}" -gt "70" ]; then
        bwa mem -M -t !{task.cpus} !{bwa_mem_opts} multifasta.fa !{fq} > bwa.sam 2> bwa-mem.stderr.txt
    else
        if [ "!{single_end}" == "true" ]; then
            bwa aln -f bwa.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[0]} > bwa-aln.stdout.txt 2> bwa-aln.stderr.txt
            bwa samse -n !{params.bwa_n} !{bwa_samse_opts} multifasta.fa bwa.sai !{fq[0]} > bwa.sam 2> bwa-samse.stderr.txt
        else
            bwa aln -f r1.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[0]} > bwa-aln.stdout.txt 2> bwa-aln.stderr.txt
            bwa aln -f r2.sai -t !{task.cpus} !{bwa_aln_opts} multifasta.fa !{fq[1]} >> bwa-aln.stdout.txt 2>> bwa-aln.stderr.txt
            bwa sampe -n !{params.bwa_n} !{bwa_sampe_opts} multifasta.fa r1.sai r2.sai !{fq[0]} !{fq[1]} > bwa.sam 2> bwa-sampe.stderr.txt
        fi
    fi

    # Write per-base coverage
    samtools view -bS bwa.sam | samtools sort -o cov.bam - > samtools.stdout.txt 2> samtools.stderr.txt
    genomeCoverageBed -ibam cov.bam -d > cov.txt 2> genomeCoverageBed.stderr.txt
    split-coverages.py mapping.txt cov.txt --outdir results

    if [[ !{params.skip_compression} == "false" ]]; then
        pigz --best -n -p !{task.cpus} results/*.txt
    fi

    # Capture version
    bwa 2>&1 | grep "Version" > bwa.version.txt 2>&1
    samtools 2>&1 | grep "Version" > samtools.version.txt 2>&1
    bedtools --version > bedtools.version.txt 2>&1
    '''

    stub:
    """
    mkdir ${PROCESS_NAME}
    mkdir mapping
    touch ${PROCESS_NAME}/${sample}
    touch mapping/${sample}
    """
}
