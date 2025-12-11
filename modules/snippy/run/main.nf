/**
 * Rapid haploid variant calling.
 *
 * This process executes snippy_run to perform analysis
 *
 * @status stable
 * @keywords variant, fastq, bacteria
 * @tags complexity:complex input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic
 * @citation snippy_run
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: List of input FastQ files of size 1 and 2 for single-end and paired-end data,
 * respectively.
 * 
 *
 * @input tuple(meta, reference)
 * - `meta`: Groovy Map containing sample information
 * - `reference`: Input file
 *
 * @output aligned_fa               A version of the reference but with - at position with depth=0 and N for 0 < depth < --mincov (does not have variants)
 * @output vcf                      The final annotated variants in VCF format
 * @output aligned_fa_error         Aligned Fa Error
 * @output vcf_error                Vcf Error
 * @output error                    Error
 * @output annotated_vcf            Annotated Vcf
 * @output bam                      The alignments in BAM format. Includes unmapped, multimapping reads. Excludes duplicates.
 * @output bai                      Index for the .bam file
 * @output bed                      The variants in BED format
 * @output consensus_fa             A version of the reference genome with all variants instantiated
 * @output consensus_subs_fa        A version of the reference genome with only substitution variants instantiated
 * @output consensus_subs_masked_fa Consensus Subs Masked Fa
 * @output coverage                 Coverage
 * @output csv                      A comma-separated version of the .tab file
 * @output filt_vcf                 The filtered variant calls from Freebayes
 * @output gff                      The variants in GFF3 format
 * @output html                     A HTML version of the .tab file
 * @output raw_vcf                  The unfiltered variant calls from Freebayes
 * @output subs_vcf                 Subs Vcf
 * @output tab                      A simple tab-separated summary of all the variants
 * @output txt                      Tab-separated columnar list of statistics
 * @output logs                     Optional tool execution logs
 * @output nf_logs                  Nextflow execution logs
 * @output versions                 Software version information (YAML format)
 */
nextflow.preview.types = true

process SNIPPY_RUN {
    tag "${prefix} - ${reference_name}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads)         : Tuple<Map, List<Path>>
    (_ref_meta, reference) : Tuple<Map, Path>

    output:
    aligned_fa               = tuple(meta, file("${prefix}.aligned.fa.gz", optional: true))
    vcf                      = tuple(meta, file("${prefix}.vcf.gz", optional: true))
    aligned_fa_error         = tuple(meta, file("${prefix}.error.aligned.fa.gz", optional: true))
    vcf_error                = tuple(meta, file("${prefix}.error.vcf.gz", optional: true))
    error                    = tuple(meta, file("${prefix}.error.txt", optional: true))
    annotated_vcf            = tuple(meta, file("${prefix}.annotated.vcf.gz"))
    bam                      = tuple(meta, file("${prefix}.bam", optional: true))
    bai                      = tuple(meta, file("${prefix}.bam.bai", optional: true))
    bed                      = tuple(meta, file("${prefix}.bed.gz"))
    consensus_fa             = tuple(meta, file("${prefix}.consensus.fa.gz"))
    consensus_subs_fa        = tuple(meta, file("${prefix}.consensus.subs.fa.gz"))
    consensus_subs_masked_fa = tuple(meta, file("${prefix}.consensus.subs.masked.fa.gz"))
    coverage                 = tuple(meta, file("${prefix}.coverage.txt.gz"))
    csv                      = tuple(meta, file("${prefix}.csv.gz"))
    filt_vcf                 = tuple(meta, file("${prefix}.filt.vcf.gz"))
    gff                      = tuple(meta, file("${prefix}.gff.gz"))
    html                     = tuple(meta, file("${prefix}.html"))
    raw_vcf                  = tuple(meta, file("${prefix}.raw.vcf.gz"))
    subs_vcf                 = tuple(meta, file("${prefix}.subs.vcf.gz"))
    tab                      = tuple(meta, file("${prefix}.tab"))
    txt                      = tuple(meta, file("${prefix}.txt"))
    logs                     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs                  = tuple(meta, files(".command.*"))
    versions                 = tuple(meta, file("versions.yml"))

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
