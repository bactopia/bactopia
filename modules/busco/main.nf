/**
 * Assess genome assembly completeness using single-copy orthologs.
 *
 * Uses [BUSCO](https://gitlab.com/ezlab/busco) (Benchmarking Universal Single-Copy Orthologs)
 * to measure the completeness of genome assemblies, gene sets, or transcriptomes by matching
 * them against a lineage-specific set of conserved orthologs.
 *
 * @status stable
 * @keywords quality control, completeness, genome, assembly, orthologs, busco
 * @tags complexity:moderate input-type:single output-type:multiple features:internet-access,resource-download
 * @citation busco
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, supplemental, results, logs, nf_logs, versions)
 * - `tsv`: Text summary report of the completeness score (C/S/D/F/M%)
 * - `supplemental`: Directory containing full tables, missing gene lists, and lineage data
 */
nextflow.preview.types = true

process BUSCO {
    tag "${prefix} - ${task.ext.busco_lineage}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, fna: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}-summary.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}-summary.txt"),
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

    def which_cat = fna.getName().endsWith(".gz") ? "zcat" : "cat"
    def fna_name = fna.getName().replace(".gz", "")
    """
    # Have to put FASTA in a directory to force batch mode in busco
    mkdir tmp-fasta
    ${which_cat} ${fna} > tmp-fasta/${fna_name}

    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    busco \\
        --cpu ${task.cpus} \\
        --in tmp-fasta/ \\
        --out supplemental \\
        --lineage ${task.ext.busco_lineage} \\
        --mode genome \\
        --download_base_url=https://busco-data2.ezlab.org/v5/data \\
        ${task.ext.args}

    # Cleanup
    find supplemental/ -name "*.log" | xargs -I {} mv {} ./
    find supplemental/ -type d -path "*logs" | xargs -I {} rm -rf {}
    find supplemental/ -type f -name "*.fna" | xargs -I {} gzip {}
    find supplemental/ -type f -name "*.faa" | xargs -I {} gzip {}
    find supplemental/ -type f -path "*hmmer_output*" -name "*.out" | xargs -I {} gzip {}
    mv supplemental/batch_summary.txt supplemental/${prefix}-summary.txt
    mv supplemental/${fna_name}/* supplemental/
    rm -rf supplemental/${fna_name} busco_downloads/ tmp*/

    # Busco outputs additional trailing tabs, clean them up
    sed -i 's/\t\t\t\$//' supplemental/${prefix}-summary.txt
    mv supplemental/${prefix}-summary.txt ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
