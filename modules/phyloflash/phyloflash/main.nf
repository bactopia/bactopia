/**
 * Rapidly reconstruct SSU rRNAs and determine phylogenetic composition.
 *
 * Uses [phyloFlash](https://github.com/HRGV/phyloFlash) to map Illumina reads to a reference
 * SSU rRNA database (typically SILVA). It then performs read-based assembly of the SSU rRNA genes
 * to accurately determine the phylogenetic composition and relative abundance of taxa in the sample.
 *
 * @status stable
 * @keywords metagenomics, taxonomy, phylogeny, ssu rrna, silva, ribosomal rna, abundance
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent
 * @citation phyloflash
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end or single-end reads in FASTQ format
 *
 * @input _silva_db
 * Path to the SILVA SSU rRNA database used for read mapping
 *
 * @input _univec_db
 * Path to the UniVec database for screening contaminants (typically optional)
 *
 * @output record(meta, supplemental, aln, summary, results, logs, nf_logs, versions)
 * - `supplemental`: Directory containing various intermediate files and detailed outputs
 * - `aln`: Reconstructed SSU rRNA sequences, aligned in FASTA format
 * - `summary`: JSON summary file containing taxonomic hits and abundance statistics
 */
nextflow.preview.types = true

process PHYLOFLASH {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, reads: Set<Path>): Record
    _silva_db       : Path
    _univec_db      : Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        supplemental: files("${prefix}/*"),
        aln: files("${prefix}/${prefix}.toalign.fasta", optional: true),
        summary: files("${prefix}/${prefix}.phyloFlash.json", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}/${prefix}.toalign.fasta", optional: true),
            files("${prefix}/${prefix}.phyloFlash.json", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    def read_opts = meta.single_end ? "-read1 ${reads.toList()[0]}" : "-read1 ${reads.toList()[0]} -read2 ${reads.toList()[1]}"
    """
    mkdir ${prefix}
    phyloFlash.pl \\
        ${task.ext.args} \\
        ${read_opts} \\
        -lib ${prefix} \\
        -dbhome . \\
        -CPUs ${task.cpus}

    jsonify-phyloflash.py ${prefix}.phyloFlash > ${prefix}.phyloFlash.json
    mv ${prefix}.* ${prefix}


    if phyloflash-summary.py ${prefix}/ | grep -q -c "WARNING: Multiple SSUs were assembled by SPAdes"; then
        MULTI="1"
    fi

    if [ "${task.ext.allow_multiple_16s}" == "true" ]; then
        MULTI="0"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyloFlash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
    END_VERSIONS
    """
}
