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
 * @citation snippy_core
 *
 * @input tuple(meta, vcf, aligned_fa)
 * - `meta`: Groovy Map containing sample information
 * - `vcf`: List of VCF files from Snippy
 * - `aligned_fa`: List of aligned FASTA files from Snippy
 *
 * @input tuple(meta, reference)
 * - `meta`: Groovy Map containing reference information
 * - `reference`: Reference genome (FASTA or GenBank format)
 *
 * @input mask
 * Optional BED file of regions to mask in the alignment
 *
 * @output supplemental   Supplemental files including individual sample alignments
 * @output aln            A core SNP alignment in FASTA format
 * @output full_aln       A whole genome SNP alignment (includes invariant sites)
 * @output clean_full_aln A whole genome SNP alignment (includes invariant sites) with Ns
 * @output tab            Tab-separated list of core SNP sites with alleles (no annotations)
 * @output vcf            Multi-sample VCF file with genotype GT tags for all discovered alleles
 * @output txt            Tab-separated list of alignment and core-size statistics
 * @output samples        List of samples included in the core alignment
 * @output logs           Optional software execution logs containing warnings/errors
 * @output nf_logs        Nextflow execution scripts and logs for debugging
 * @output versions       A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SNIPPY_CORE {
    tag "${prefix}"
    label "process_medium"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, _vcf, _aligned_fa) : Tuple<Map, Path, Path>
    (_ref_meta, reference)   : Tuple<Map, Path>
    mask                     : List<Path>

    output:
    supplemental   = tuple(meta, files("snippy-core/*"))
    aln            = tuple(meta, file("snippy-core/${prefix}.aln.gz"))
    full_aln       = tuple(meta, file("${prefix}.full.aln.gz"))
    clean_full_aln = tuple(meta, file("${prefix}-clean.full.aln.gz"))
    tab            = tuple(meta, file("snippy-core/${prefix}.tab.gz"))
    vcf            = tuple(meta, file("snippy-core/${prefix}.vcf.gz"))
    txt            = tuple(meta, file("snippy-core/${prefix}.txt"))
    samples        = tuple(meta, file("${reference_name}.samples.txt"))
    logs           = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs        = tuple(meta, files(".command.*"))
    versions       = tuple(meta, file("versions.yml"))

    script:
    reference_name = reference.getSimpleName()
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = ""
    meta.logs_dir = "${meta.process_name}/logs"
    def mask_opt = mask.size() == 1 ? "--mask ${mask[0].getName()}" : ""
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

    # Cleanup the alignment
    snippy-clean_full_aln ${prefix}.full.aln > ${prefix}-clean.full.aln
    rm -rf *ref.fa ${final_reference} samples/

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
