/**
 * Taxonomic classification of bacterial and archaeal genomes using GTDB-Tk.
 *
 * Uses [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) to assign objective taxonomic
 * classifications to genome assemblies based on the [Genome Taxonomy Database](https://gtdb.ecogenomic.org/).
 * It identifies marker genes, aligns them, and places the genome into the reference tree to determine taxonomy.
 *
 * @status stable
 * @keywords taxonomy, classification, phylogeny, gtdb, bacteria, archaea, marker genes
 * @tags complexity:complex input-type:multiple output-type:multiple features:database-dependent,conditional-logic
 * @citation gtdb_tk
 *
 * @note Database Required
 * Requires the massive GTDB-Tk database (~60GB+) to be available.
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input db
 * Path (or Set of paths) to the GTDB-Tk reference database
 *
 * @output record(meta, bac_tsv, ar_tsv, results, logs, nf_logs, versions)
 * - `bac_tsv`: The bacterial classification summary file containing the taxonomic assignment
 * - `ar_tsv`: The archaeal classification summary file containing the taxonomic assignment
 *
 * @results supplemental
 * - `align/`: Multiple sequence alignments of identified marker genes
 * - `classify/`: Detailed classification results and tree placement files
 * - `identify/`: Marker gene identification results and translation tables
 */
nextflow.preview.types = true

process GTDBTK_CLASSIFYWF {
    tag "${prefix}"
    label 'process_high'
    label 'process_high_memory'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, fna: Set<Path>): Record
    db: Path

    stage:
    stageAs 'fna-tmp/*', fna
    stageAs 'gtdb/*', db

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        bac_tsv: file("${prefix}.bac120.summary.tsv", optional: true),
        ar_tsv: file("${prefix}.ar53.summary.tsv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.bac120.summary.tsv", optional: true),
            files("${prefix}.ar53.summary.tsv", optional: true),
            files("supplemental/*")
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

    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        export GTDBTK_DATA_PATH="\$(realpath \$(find database/ -path "*metadata*" -name "metadata.txt" | sed 's=/metadata/metadata.txt=='))"
    else
        export GTDBTK_DATA_PATH="\$(readlink ${db})"
    fi
    mkdir fna
    cp -L fna-tmp/* fna/
    find fna/ -name "*.fna.gz" | xargs gunzip

    gtdbtk classify_wf \\
        ${task.ext.args} \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${task.cpus} \\
        --genome_dir ./fna \\
        --out_dir supplemental \\
        --skip_ani_screen \\
        --prefix ${prefix}
    mv supplemental/*.log ./
    mv supplemental/classify/*.summary.tsv ./

    # Cleanup
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database/
    fi
    if [ "${task.ext.gtdb_keep_msa}" == "false" ]; then
        rm -rf supplemental/align/*.msa.fasta.gz
    fi
    rm -rf fna/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb-tk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
