/**
 * Rapid haploid variant calling and core genome alignment.
 *
 * Uses [Snippy](https://github.com/tseemann/snippy) to find SNPs and indels between a haploid
 * reference genome and your Next-Generation Sequencing (NGS) sequence reads. It maps reads to
 * the reference, calls variants, and generates a consensus sequence.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se) where each read slot is Path?
 *
 * @status stable
 * @keywords snippy, variant calling, snp, indel, alignment, bacteria
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation snippy
 *
 * @input record(meta, r1, r2, se)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 *
 * @input record(meta, reference)
 * - `meta`: Groovy Map containing reference information
 * - `reference`: Reference genome (FASTA or GenBank format)
 *
 * @output record(meta, aligned_fa, vcf, aligned_fa_error, vcf_error, error, annotated_vcf, bam, bai, bed, consensus_fa, consensus_subs_fa, consensus_subs_masked_fa, coverage, csv, filt_vcf, gff, html, raw_vcf, subs_vcf, tab, txt, results, logs, nf_logs, versions)
 * - `aligned_fa`: A version of the reference with - at zero coverage positions
 * - `vcf`: The final annotated variants in VCF format
 * - `aligned_fa_error`: Aligned FASTA file generated during error state
 * - `vcf_error`: VCF file generated during error state
 * - `error`: Error log text file
 * - `annotated_vcf`: Annotated VCF file
 * - `bam`: The alignments in BAM format (includes unmapped/multimapping)
 * - `bai`: Index for the BAM file
 * - `bed`: The variants in BED format
 * - `consensus_fa`: Reference genome with all variants instantiated
 * - `consensus_subs_fa`: Reference genome with only substitution variants instantiated
 * - `consensus_subs_masked_fa`: Reference genome with substitutions instantiated and low coverage masked
 * - `coverage`: Per-base coverage depth information
 * - `csv`: A comma-separated summary of variants
 * - `filt_vcf`: The filtered variant calls from Freebayes
 * - `gff`: The variants in GFF3 format
 * - `html`: A HTML summary of the variants
 * - `raw_vcf`: The unfiltered variant calls from Freebayes
 * - `subs_vcf`: VCF containing only substitution variants
 * - `tab`: A simple tab-separated summary of all variants
 * - `txt`: Tab-separated columnar list of alignment statistics
 */
nextflow.preview.types = true

process SNIPPY_RUN {
    tag "${prefix} - ${reference_name}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?): Record
    reference: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        aligned_fa: file("${prefix}.aligned.fa.gz", optional: true),
        vcf: file("${prefix}.vcf.gz", optional: true),
        aligned_fa_error: file("${prefix}.error.aligned.fa.gz", optional: true),
        vcf_error: file("${prefix}.error.vcf.gz", optional: true),
        error: file("${prefix}.error.txt", optional: true),
        annotated_vcf: file("${prefix}.annotated.vcf.gz"),
        bam: file("${prefix}.bam", optional: true),
        bai: file("${prefix}.bam.bai", optional: true),
        bed: file("${prefix}.bed.gz"),
        consensus_fa: file("${prefix}.consensus.fa.gz"),
        consensus_subs_fa: file("${prefix}.consensus.subs.fa.gz"),
        consensus_subs_masked_fa: file("${prefix}.consensus.subs.masked.fa.gz"),
        coverage: file("${prefix}.coverage.txt.gz"),
        csv: file("${prefix}.csv.gz"),
        filt_vcf: file("${prefix}.filt.vcf.gz"),
        gff: file("${prefix}.gff.gz"),
        html: file("${prefix}.html"),
        raw_vcf: file("${prefix}.raw.vcf.gz"),
        subs_vcf: file("${prefix}.subs.vcf.gz"),
        tab: file("${prefix}.tab"),
        txt: file("${prefix}.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.aligned.fa.gz", optional: true),
            files("${prefix}.vcf.gz", optional: true),
            files("${prefix}.error.aligned.fa.gz", optional: true),
            files("${prefix}.error.vcf.gz", optional: true),
            files("${prefix}.error.txt", optional: true),
            files("${prefix}.annotated.vcf.gz"),
            files("${prefix}.bam", optional: true),
            files("${prefix}.bam.bai", optional: true),
            files("${prefix}.bed.gz"),
            files("${prefix}.consensus.fa.gz"),
            files("${prefix}.consensus.subs.fa.gz"),
            files("${prefix}.consensus.subs.masked.fa.gz"),
            files("${prefix}.coverage.txt.gz"),
            files("${prefix}.csv.gz"),
            files("${prefix}.filt.vcf.gz"),
            files("${prefix}.gff.gz"),
            files("${prefix}.html"),
            files("${prefix}.raw.vcf.gz"),
            files("${prefix}.subs.vcf.gz"),
            files("${prefix}.tab"),
            files("${prefix}.txt")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    meta.single_end = has_se && !has_r1 && !has_r2

    // Build read inputs for snippy
    def read_inputs = meta.single_end ? "--se ${se}" : "--R1 ${r1} --R2 ${r2}"
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

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${final_reference}
    fi
    rm -rf tmp_snippy/ ${prefix}/reference ${prefix}/ref.fa* ${prefix}/${prefix}.vcf.gz*

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
    rm -rf ${prefix}/

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
