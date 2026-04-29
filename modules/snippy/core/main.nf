/**
 * Core-SNP alignment from Snippy outputs.
 *
 * Uses [Snippy](https://github.com/tseemann/snippy) to generate a core genome alignment
 * from multiple Snippy outputs. It combines variant calls (VCF) and alignments to produce
 * a core SNP alignment, which can be used for phylogenetic analysis.
 *
 * @status stable
 * @keywords snippy, core genome, alignment, phylogeny, snp, bacteria
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation snippy, bcftools, bedtools, freebayes, seqtk, snpeff, vcf_annotator, vcflib, vt
 *
 * @input record(meta, _vcf, _aligned_fa)
 * - `meta`: Groovy Record containing sample information
 * - `_vcf`: List of VCF files from Snippy
 * - `_aligned_fa`: List of aligned FASTA files from Snippy
 *
 * @input record(meta, reference)
 * - `meta`: Groovy Record containing reference information
 * - `reference`: Reference genome (FASTA or GenBank format)
 *
 * @input mask?
 * BED file of regions to mask in the alignment
 *
 * @output record(meta, supplemental, aln, full_aln, clean_full_aln, tab, vcf, txt, samples, results, logs, nf_logs, versions)
 * - `supplemental`: Supplemental files including individual sample alignments
 * - `aln`: A core SNP alignment in FASTA format
 * - `full_aln`: A whole genome SNP alignment (includes invariant sites)
 * - `clean_full_aln`: A whole genome SNP alignment (includes invariant sites) with Ns
 * - `tab`: Tab-separated list of core SNP sites with alleles (no annotations)
 * - `vcf`: Multi-sample VCF file with genotype GT tags for all discovered alleles
 * - `txt`: Tab-separated list of alignment and core-size statistics
 * - `samples`: List of samples included in the core alignment
 */
nextflow.enable.types = true

process SNIPPY_CORE {
    tag "${prefix}"
    label "process_medium"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        _vcf: Set<Path>,
        _aligned_fa: Set<Path>
    )
    reference: Path
    mask: Path?

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        supplemental: files("snippy-core/*"),
        aln: file("snippy-core/${prefix}.aln.gz"),
        full_aln: file("${prefix}.full.aln.gz"),
        clean_full_aln: file("${prefix}-clean.full.aln.gz"),
        tab: file("snippy-core/${prefix}.tab.gz"),
        vcf: file("snippy-core/${prefix}.vcf.gz"),
        txt: file("snippy-core/${prefix}.txt"),
        samples: file("${reference_name}.samples.txt"),
        // Generic fields (used for publishing)
        results: [
            files("snippy-core/*"),
            files("snippy-core/${prefix}.tab.gz"),
            files("snippy-core/${prefix}.vcf.gz"),
            files("snippy-core/${prefix}.txt"),
            files("snippy-core/${prefix}.aln.gz"),
            files("${prefix}.full.aln.gz"),
            files("${prefix}-clean.full.aln.gz"),
            files("${reference_name}.samples.txt")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    reference_name = reference.getSimpleName()
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        process_name: task.ext.process_name,
        output_dir: "",
        logs_dir: "${task.ext.process_name}/logs"
    )

    // Tool specific variables
    def mask_opt = mask != null ? "--mask ${mask.getName()}" : ""
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    def final_reference = reference.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${reference} > ${final_reference}
    fi

    # Collect samples into necessary folders
    mkdir samples
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//;s=./==' > ${reference_name}.samples.txt
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//' | xargs -I {} bash -c 'mkdir samples/{}'
    find . -name "*.vcf.gz" | sed 's/.vcf.gz//' | xargs -I {} bash -c 'gzip -cdf {}.vcf.gz > samples/{}/{}.vcf'
    find . -name "*.aligned.fa.gz" | sed 's/.aligned.fa.gz//' | xargs -I {} bash -c 'gzip -cdf {}.aligned.fa.gz > samples/{}/{}.aligned.fa'

    # Run snippy-core
    snippy-core \\
        ${task.ext.args} \\
        --ref ${final_reference} \\
        --prefix ${prefix} \\
        ${mask_opt} \\
        samples/*

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${final_reference}
    fi
    snippy-clean_full_aln ${prefix}.full.aln > ${prefix}-clean.full.aln
    rm -rf *ref.fa samples/

    # Compress outputs
    if [[ ${task.ext.skip_compression} == "false" ]]; then
        pigz -n --best -p ${task.cpus} ${prefix}.aln
        pigz -n --best -p ${task.cpus} ${prefix}.full.aln
        pigz -n --best -p ${task.cpus} ${prefix}-clean.full.aln
        pigz -n --best -p ${task.cpus} ${prefix}.tab
        pigz -n --best -p ${task.cpus} ${prefix}.vcf
    fi

    # Move outputs
    mkdir snippy-core
    mv ${prefix}* snippy-core/
    mv snippy-core/${prefix}-clean.full.aln.gz snippy-core/${prefix}.full.aln.gz ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy //')
    END_VERSIONS
    """
}
