/**
 * Rapid haploid variant calling and core genome alignment.
 *
 * Uses [Snippy](https://github.com/tseemann/snippy) to find SNPs and indels between a haploid
 * reference genome and your Next-Generation Sequencing (NGS) sequence reads. It maps reads to
 * the reference, calls variants, and generates a consensus sequence.
 *
 * @status stable
 * @keywords snippy, variant calling, snp, indel, alignment, bacteria
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation snippy
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: FASTQ reads (single-end or paired-end)
 *
 * @input tuple(meta, reference)
 * - `meta`: Groovy Map containing reference information
 * - `reference`: Reference genome (FASTA or GenBank format)
 *
 * @output aligned_fa               A version of the reference with - at zero coverage positions
 * @output vcf                      The final annotated variants in VCF format
 * @output aligned_fa_error         Aligned FASTA file generated during error state
 * @output vcf_error                VCF file generated during error state
 * @output error                    Error log text file
 * @output annotated_vcf            Annotated VCF file
 * @output bam                      The alignments in BAM format (includes unmapped/multimapping)
 * @output bai                      Index for the BAM file
 * @output bed                      The variants in BED format
 * @output consensus_fa             Reference genome with all variants instantiated
 * @output consensus_subs_fa        Reference genome with only substitution variants instantiated
 * @output consensus_subs_masked_fa Reference genome with substitutions instantiated and low coverage masked
 * @output coverage                 Per-base coverage depth information
 * @output csv                      A comma-separated summary of variants
 * @output filt_vcf                 The filtered variant calls from Freebayes
 * @output gff                      The variants in GFF3 format
 * @output html                     A HTML summary of the variants
 * @output raw_vcf                  The unfiltered variant calls from Freebayes
 * @output subs_vcf                 VCF containing only substitution variants
 * @output tab                      A simple tab-separated summary of all variants
 * @output txt                      Tab-separated columnar list of alignment statistics
 * @output logs                     Optional software execution logs containing warnings/errors
 * @output nf_logs                  Nextflow execution scripts and logs for debugging
 * @output versions                 A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SNIPPY_RUN {
    tag "${prefix} - ${reference_name}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads)         : Tuple<Map, Set<Path>>
    (_ref_meta, reference) : Tuple<Map, Set<Path>>

    output:
    aligned_fa               = tuple(meta, files("${prefix}.aligned.fa.gz", optional: true))
    vcf                      = tuple(meta, files("${prefix}.vcf.gz", optional: true))
    aligned_fa_error         = tuple(meta, files("${prefix}.error.aligned.fa.gz", optional: true))
    vcf_error                = tuple(meta, files("${prefix}.error.vcf.gz", optional: true))
    error                    = tuple(meta, files("${prefix}.error.txt", optional: true))
    annotated_vcf            = tuple(meta, files("${prefix}.annotated.vcf.gz"))
    bam                      = tuple(meta, files("${prefix}.bam", optional: true))
    bai                      = tuple(meta, files("${prefix}.bam.bai", optional: true))
    bed                      = tuple(meta, files("${prefix}.bed.gz"))
    consensus_fa             = tuple(meta, files("${prefix}.consensus.fa.gz"))
    consensus_subs_fa        = tuple(meta, files("${prefix}.consensus.subs.fa.gz"))
    consensus_subs_masked_fa = tuple(meta, files("${prefix}.consensus.subs.masked.fa.gz"))
    coverage                 = tuple(meta, files("${prefix}.coverage.txt.gz"))
    csv                      = tuple(meta, files("${prefix}.csv.gz"))
    filt_vcf                 = tuple(meta, files("${prefix}.filt.vcf.gz"))
    gff                      = tuple(meta, files("${prefix}.gff.gz"))
    html                     = tuple(meta, files("${prefix}.html"))
    raw_vcf                  = tuple(meta, files("${prefix}.raw.vcf.gz"))
    subs_vcf                 = tuple(meta, files("${prefix}.subs.vcf.gz"))
    tab                      = tuple(meta, files("${prefix}.tab"))
    txt                      = tuple(meta, files("${prefix}.txt"))
    logs                     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs                  = tuple(meta, files(".command.*"))
    versions                 = tuple(meta, files("versions.yml"))

    script:
    reference_name = reference.toList()[0].getSimpleName()
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${reference_name}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${reference_name}/logs"
    meta.process_name = task.ext.process_name
    def read_inputs = _meta.single_end ? "--se ${reads.toList()[0]}" : "--R1 ${reads.toList()[0]} --R2 ${reads.toList()[1]}"
    def is_compressed = reference.toList()[0].getName().endsWith(".gz") ? true : false
    def final_reference = reference.toList()[0].getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${reference} > ${final_reference}
    fi

    if ! head -n 1 ${final_reference} | grep "^LOCUS"; then
        echo "ERROR: Reference file (${reference}) does not appear to be a GenBank file"
        exit 1
    fi

    mkdir tmp_snippy/
    snippy \\
        ${task.ext.args} \\
        --cpus ${task.cpus} \\
        --outdir ${prefix} \\
        --reference ${final_reference} \\
        --prefix ${prefix} \\
        --tmpdir tmp_snippy/ \\
        ${read_inputs}

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
