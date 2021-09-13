nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "call_variants"

process CALL_VARIANTS {
    /*
    Identify variants (SNPs/InDels) against a set of reference genomes
    using Snippy.

    If applicable, download the nearest RefSeq genomes (based on Mash) to have variants called against.

    Exitcode 75 is due to being unable to download from NCBI (e.g. FTP down at the time)
    Downloads will be attempted 300 times total before giving up. On failure to download
    variants will not be called against the nearest completed genome.
    */
    tag "${sample} - ${reference_name}"
    label "max_cpu_75"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME, logs_subdir:reference_name) }

    input:
    tuple val(sample), val(single_end), path(fq)
    each path(reference)

    output:
    path "results/*"
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    PROCESS_NAME = "call_variants"
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getSimpleName()
    fastq = single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
    bwaopt = params.bwaopt ? "--bwaopt 'params.bwaopt'" : ""
    fbopt = params.fbopt ? "--fbopt 'params.fbopt'" : ""
    no_cache = params.no_cache ? '-N' : ''
    tie_break = params.random_tie_break ? "--random_tie_break" : ""
    '''
    REFERENCE="!{reference}"
    REFERENCE_NAME="!{reference_name}"
    if [[ "!{reference}" == "refseq-genomes.msh" ]]; then
        # Get Mash distance (For paired-end, use just R1)
        mash dist -k 31 -s 10000 -t !{fq[0]} ${REFERENCE} | grep -v "query" | sort -k 2,2 > distances.txt

        # Pick top genome and download
        printf "accession\\tdistance\\tlatest_accession\\tupdated\\n" > mash-dist.txt
        select-references.py distances.txt 1 !{tie_break} >> mash-dist.txt
        grep -v distance mash-dist.txt | cut -f3 > download-list.txt
        ncbi-genome-download bacteria -l complete -o ./ -F genbank -p !{task.cpus} -A download-list.txt \
            -r !{params.max_retry} !{no_cache} > ncbi-genome-download.stdout.txt 2> ncbi-genome-download.stderr.txt

        # Move and uncompress genomes
        mkdir genbank_temp
        find refseq -name "*.gbff.gz" | xargs -I {} mv {} genbank_temp/
        rename 's/(GC[AF]_\\d+).*/$1/' genbank_temp/*
        ls genbank_temp/ | xargs -I {} sh -c 'gzip -cd genbank_temp/{} > {}.gbk'
        rm -rf genbank_temp/ refseq/

        # Update variables
        REFERENCE=$(ls *.gbk)
        REFERENCE_NAME=$(basename ${REFERENCE} .gbk)

        # Capture version
        mash --version > mash.version.txt 2>&1
        ncbi-genome-download --version > ncbi-genome-download.version.txt 2>&1
    fi

    snippy !{fastq} \
        --ref ${REFERENCE} \
        --cpus !{task.cpus} \
        --ram !{snippy_ram} \
        --outdir ${REFERENCE_NAME} \
        --prefix !{sample} \
        --mapqual !{params.mapqual} \
        --basequal !{params.basequal} \
        --mincov !{params.mincov} \
        --minfrac !{params.minfrac} \
        --minqual !{params.minqual} \
        --maxsoft !{params.maxsoft} !{bwaopt} !{fbopt} > snippy.stdout.txt 2> snippy.stderr.txt

    # Add GenBank annotations
    vcf-annotator ${REFERENCE_NAME}/!{sample}.vcf ${REFERENCE} > ${REFERENCE_NAME}/!{sample}.annotated.vcf 2> vcf-annotator.stderr.txt

    # Get per-base coverage
    grep "^##contig" ${REFERENCE_NAME}/!{sample}.vcf > ${REFERENCE_NAME}/!{sample}.full-coverage.txt
    genomeCoverageBed -ibam ${REFERENCE_NAME}/!{sample}.bam -d >> ${REFERENCE_NAME}/!{sample}.full-coverage.txt 2> genomeCoverageBed.stderr.txt
    cleanup-coverage.py ${REFERENCE_NAME}/!{sample}.full-coverage.txt > ${REFERENCE_NAME}/!{sample}.coverage.txt
    rm ${REFERENCE_NAME}/!{sample}.full-coverage.txt

    # Mask low coverage regions
    mask-consensus.py !{sample} ${REFERENCE_NAME} \
        ${REFERENCE_NAME}/!{sample}.consensus.subs.fa \
        ${REFERENCE_NAME}/!{sample}.subs.vcf \
        ${REFERENCE_NAME}/!{sample}.coverage.txt \
        --mincov !{params.mincov} > ${REFERENCE_NAME}/!{sample}.consensus.subs.masked.fa 2> mask-consensus.stderr.txt

    # Clean Up
    rm -rf ${REFERENCE_NAME}/reference ${REFERENCE_NAME}/ref.fa* ${REFERENCE_NAME}/!{sample}.vcf.gz*

    if [[ "!{reference}" == "refseq-genomes.msh" ]]; then
        mv distances.txt ${REFERENCE_NAME}/mash-distances.txt
        mv ${REFERENCE} ${REFERENCE_NAME}/
    fi

    if [[ !{params.skip_compression} == "false" ]]; then
        find ${REFERENCE_NAME}/ -type f | \
            grep -v -E "\\.bam$|\\.bai$|\\.log$|\\.txt$|\\.html$|\\.tab$" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
        pigz -n --best -p !{task.cpus} ${REFERENCE_NAME}/!{sample}.coverage.txt
    fi

    mkdir results
    mv ${REFERENCE_NAME}/ results/

    # Capture verisons
    snippy --version > snippy.version.txt 2>&1
    vcf-annotator --version > vcf-annotator.version.txt 2>&1
    bedtools --version > bedtools.version.txt 2>&1
    '''

    stub:
    reference_name = reference.getSimpleName()
    """
    mkdir ${reference_name}
    touch ${reference_name}/*
    """
}
