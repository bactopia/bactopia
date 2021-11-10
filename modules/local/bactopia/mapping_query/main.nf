nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
params.options = [:]
options = initOptions(params.options, 'mapping_query')

process MAPPING_QUERY {
    /*
    Map FASTQ reads against a given set of FASTA files using BWA.
    */
    tag "${meta.id}"
    label "max_cpus"
    label "mapping_query"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    input:
    tuple val(meta), path(fq)
    each path(query)

    output:
    path "results/*"
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

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
        if [ "!{meta.single_end}" == "true" ]; then
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

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    mapping_query:
        bedtools: $(echo $(bedtools --version 2>&1) | sed 's/bedtools v//')
        bwa: $(echo $(bwa 2>&1) | sed 's/^.*Version: //;s/ .*$//')
        pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
        samtools: $(echo $(samtools --version 2>&1) |sed 's/^.*samtools //;s/ .*$//')
    END_VERSIONS
    '''
}
