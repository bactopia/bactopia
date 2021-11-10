nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
params.options = [:]
options = initOptions(params.options, 'call_variants')

process CALL_VARIANTS {
    /*
    Identify variants (SNPs/InDels) against a set of reference genomes
    using Snippy.

    If applicable, download the nearest RefSeq genomes (based on Mash) to have variants called against.

    Exitcode 75 is due to being unable to download from NCBI (e.g. FTP down at the time)
    Downloads will be attempted 300 times total before giving up. On failure to download
    variants will not be called against the nearest completed genome.
    */
    tag "${meta.id} - ${reference_name}"
    label "max_cpu_75"
    label "call_variants"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options, logs_subdir:reference_name) }

    input:
    tuple val(meta), path(fq)
    each path(reference)

    output:
    path "results/*"
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
    snippy_ram = task.memory.toString().split(' ')[0]
    reference_name = reference.getSimpleName()
    fastq = meta.single_end ? "--se ${fq[0]}" : "--R1 ${fq[0]} --R2 ${fq[1]}"
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
        --prefix !{meta.id} \
        --mapqual !{params.mapqual} \
        --basequal !{params.basequal} \
        --mincov !{params.mincov} \
        --minfrac !{params.minfrac} \
        --minqual !{params.minqual} \
        --maxsoft !{params.maxsoft} !{bwaopt} !{fbopt} > snippy.stdout.txt 2> snippy.stderr.txt
    mv ${REFERENCE_NAME}/!{meta.id}.log ./

    # Add GenBank annotations
    vcf-annotator ${REFERENCE_NAME}/!{meta.id}.vcf ${REFERENCE} > ${REFERENCE_NAME}/!{meta.id}.annotated.vcf 2> vcf-annotator.stderr.txt

    # Get per-base coverage
    grep "^##contig" ${REFERENCE_NAME}/!{meta.id}.vcf > ${REFERENCE_NAME}/!{meta.id}.full-coverage.txt
    genomeCoverageBed -ibam ${REFERENCE_NAME}/!{meta.id}.bam -d >> ${REFERENCE_NAME}/!{meta.id}.full-coverage.txt 2> genomeCoverageBed.stderr.txt
    cleanup-coverage.py ${REFERENCE_NAME}/!{meta.id}.full-coverage.txt > ${REFERENCE_NAME}/!{meta.id}.coverage.txt
    rm ${REFERENCE_NAME}/!{meta.id}.full-coverage.txt

    # Mask low coverage regions
    mask-consensus.py !{meta.id} ${REFERENCE_NAME} \
        ${REFERENCE_NAME}/!{meta.id}.consensus.subs.fa \
        ${REFERENCE_NAME}/!{meta.id}.subs.vcf \
        ${REFERENCE_NAME}/!{meta.id}.coverage.txt \
        --mincov !{params.mincov} > ${REFERENCE_NAME}/!{meta.id}.consensus.subs.masked.fa 2> mask-consensus.stderr.txt

    # Clean Up
    rm -rf ${REFERENCE_NAME}/reference ${REFERENCE_NAME}/ref.fa* ${REFERENCE_NAME}/!{meta.id}.vcf.gz*

    if [[ "!{reference}" == "refseq-genomes.msh" ]]; then
        mv distances.txt ${REFERENCE_NAME}/mash-distances.txt
        mv ${REFERENCE} ${REFERENCE_NAME}/
    fi

    if [[ !{params.skip_compression} == "false" ]]; then
        find ${REFERENCE_NAME}/ -type f | \
            grep -v -E "\\.bam$|\\.bai$|\\.log$|\\.txt$|\\.html$|\\.tab$" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
        pigz -n --best -p !{task.cpus} ${REFERENCE_NAME}/!{meta.id}.coverage.txt
    fi

    mkdir results
    mv ${REFERENCE_NAME}/ results/

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    call_variants:
        bedtools: $(echo $(bedtools --version 2>&1) | sed 's/bedtools v//')
        mash: $(echo $(mash --version 2>&1))
        ncbi-genome-download: $(echo $(ncbi-genome-download --version 2>&1))
        pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
        snippy: $(echo $(snippy --version 2>&1) | sed 's/snippy //')
        vcf-annotator: $(echo $(vcf-annotator --version 2>&1) | sed 's/vcf-annotator.py //')
    END_VERSIONS
    '''
}
