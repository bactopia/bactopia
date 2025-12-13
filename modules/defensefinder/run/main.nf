/**
 * Detect anti-phage defense systems using HMM profiles.
 *
 * Uses [DefenseFinder](https://github.com/mdmparis/defense-finder) to systematically search
 * protein sequences for known antiviral defense systems (e.g., CRISPR-Cas, Restriction-Modification,
 * TA systems, CBASS) using MacSyFinder and a dedicated database of HMM models.
 *
 * @status stable
 * @keywords bacteria, defense systems, antiviral, phage, crispr, restriction-modification, hmm, macsyfinder
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,compression
 * @citation defensefinder
 *
 * @note Database Required
 * Requires the DefenseFinder HMM database to be available.
 *
 * @input tuple(meta, proteins)
 * - `meta`: Groovy Map containing sample information
 * - `proteins`: Protein sequences in FASTA format (amino acids)
 *
 * @input db
 * Directory containing the DefenseFinder models database
 *
 * @output genes_tsv     Tab-delimited list of detected defense genes
 * @output hmmer_tsv     Tab-delimited list of HMMER hits used for detection
 * @output systems_tsv   Tab-delimited summary of detected defense systems
 * @output proteins      Protein sequences of the detected defense genes (*.prt)
 * @output proteins_index Index file for the protein sequences (*.prt.idx)
 * @output macsydata_raw compressed tarball of raw MacSyFinder data (optional)
 * @output logs          Optional software execution logs containing warnings/errors
 * @output nf_logs       Nextflow execution scripts and logs for debugging
 * @output versions      A YAML formatted file with software versions
 */
nextflow.preview.types = true

process DEFENSEFINDER_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, proteins) : Tuple<Map, Set<Path>>
    db             : Path

    output:
    genes_tsv      = tuple(meta, files("*_defense_finder_genes.tsv"))
    hmmer_tsv      = tuple(meta, files("*_defense_finder_hmmer.tsv"))
    systems_tsv    = tuple(meta, files("*_defense_finder_systems.tsv"))
    proteins       = tuple(meta, files("*.prt"))
    proteins_index = tuple(meta, files("*.prt.idx"))
    macsydata_raw  = tuple(meta, files("${prefix}.macsydata.tar.gz", optional: true))
    logs           = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs        = tuple(meta, files(".command.*"))
    versions       = tuple(meta, files("versions.yml"))

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
    """
    # Extract database
    # Use custom TMPDIR to prevent FileExistsError related to writing to same tmpdir (/tmp/tmp-macsy-cache/)
    tar -xf ${db}
    mkdir -p df-tmp/df
    TMPDIR=df-tmp/df HOME=df-tmp/ macsydata \\
        install \\
        --target defense-finder/ \\
        models/defense-finder-models-v${task.ext.df_models_version}.tar.gz

    mkdir -p df-tmp/cf
    TMPDIR=df-tmp/cf HOME=df-tmp/ macsydata \\
        install \\
        --target defense-finder/ \\
        models/CasFinder-${task.ext.casfinder_version}.tar.gz

    TMPDIR=df-tmp/ HOME=df-tmp/ defense-finder \\
        run \\
        ${task.ext.args} \\
        --workers ${task.cpus} \\
        --models-dir defense-finder/ \\
        ${proteins}

    if [ "${task.ext.df_preserveraw}" == "true" ]; then
        tar -czf ${prefix}.macsydata.tar.gz defense-finder-tmp/
        rm -rf defense-finder-tmp/ 
    fi

    # Cleanup intermediate files and unused outputs
    rm -rf models/ defense-finder/ df-tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: ${task.ext.df_version}
        defense-finder-models: ${task.ext.df_models_version}
        casfinder-models: ${task.ext.casfinder_version}
    END_VERSIONS
    """
}
