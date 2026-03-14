/**
 * Create genomic sketches and perform rapid taxonomic classification.
 *
 * Uses [Mash](https://mash.readthedocs.io/) and [Sourmash](https://sourmash.readthedocs.io/) to
 * create MinHash sketches of the input sequences. These sketches are then queried against
 * pre-built databases ([RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) and
 * [GTDB](https://gtdb.ecogenomic.org/) to identify the closest reference genomes.
 *
 * @status stable
 * @keywords bacteria, taxonomy, classification, minhash, sketch, mash, sourmash, refseq, gtdb
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,compression
 * @citation mash, sourmash
 *
 * @note Databases Required
 * Requires the pre-compiled RefSeq (Mash) and GTDB (Sourmash) databases, usually downloaded
 * by the `datasets` module.
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input mash_db
 * Path to the Mash RefSeq database 
 *
 * @input sourmash_db
 * Path to the Sourmash GTDB LCA database
 *
 * @output record(meta, sig, msh, mash, sourmash, results, logs, nf_logs, versions)
 * - `sig`: The Sourmash signature file (*.sig)
 * - `msh`: The Mash sketch files for k=21 and k=31 (*.msh)
 * - `mash`: A classification report of Mash Screen results against RefSeq database
 * - `sourmash`: A classification report from Sourmash LCA against GTDB database
 */
nextflow.preview.types = true

process SKETCHER {
    tag "${prefix}"
    label "process_low"

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, assembly: Path): Record
    mash_db    : Path
    sourmash_db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        sig: file("${prefix}.sig"),
        msh: files("${prefix}-k*.msh"),
        mash: file("${prefix}-mash-refseq88-k21.txt"),
        sourmash: file("${prefix}-sourmash-gtdb-rs207-k31.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.sig"),
            files("${prefix}-k*.msh"),
            files("${prefix}-mash-refseq88-k21.txt"),
            files("${prefix}-sourmash-gtdb-rs207-k31.txt")
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
    meta.output_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/main/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = mash_db.getName().endsWith(".xz") ? true : false
    def mash_name = mash_db.getName().replace(".xz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        xz -c -d ${mash_db} > ${mash_name}
    fi

    gzip -cd ${assembly} | mash sketch -o ${prefix}-k21 -k 21 ${task.ext.args} -I ${prefix} -
    gzip -cd ${assembly} | mash sketch -o ${prefix}-k31 -k 31 ${task.ext.args} -I ${prefix} -
    sourmash sketch dna ${task.ext.args2} --merge ${prefix} -o ${prefix}.sig ${assembly}

    # Mash Screen
    echo "identity<TAB>shared-hashes<TAB>median-multiplicity<TAB>p-value<TAB>query-ID<TAB>query-comment" | sed 's/<TAB>/\t/g' > ${prefix}-mash-refseq88-k21.txt
    gzip -cd ${assembly} | mash screen ${task.ext.args3} -p ${task.cpus} ${mash_name} - | sort -gr >> ${prefix}-mash-refseq88-k21.txt

    # Sourmash classify
    sourmash lca classify --query ${prefix}.sig --db ${sourmash_db} > ${prefix}-sourmash-gtdb-rs207-k31.txt

    # Cleanup
    rm -rf ${mash_name}

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/.*sourmash //;')
    END_VERSIONS
    """
}
