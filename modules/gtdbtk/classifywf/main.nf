/**
 * GTDB-Tk classification workflow for bacterial and archaeal genomes.
 *
 * This process executes gtdbtk_classifywf to perform analysis
 *
 * @status stable
 * @keywords gtdb, taxonomy, classification, archaea, bacteria
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, database-dependent
 * @citation gtdbtk_classifywf
 *
 * @note Requires external database to be available
 *
 * @input tuple(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Genome assemblies in FASTA format
 *
 * @input db
 * GTDB database directory or tarball
 *
 * @output supplemental Supplemental
 * @output tsv          Summary classification file
 * @output logs         Optional tool execution logs
 * @output nf_logs      Nextflow execution logs
 * @output versions     Software version information (YAML format)
 */
nextflow.preview.types = true

process GTDBTK_CLASSIFYWF {
    tag "${prefix}"
    label 'process_high'
    label 'process_high_memory'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fna) : Tuple<Map, Set<Path>>
    db           : Set<Path>

    stage:
    stageAs 'fna-tmp/*', fna
    stageAs 'gtdb/*', db

    output:
    supplemental = tuple(meta, files("supplemental/*"))
    tsv          = tuple(meta, files("${prefix}.*.summary.tsv"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, file("versions.yml"))

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
    def is_tarball = db.toList()[0].getName().endsWith(".tar.gz") ? true : false
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
    mv supplemental/*.summary.tsv ./

    # Cleanup
    if [ "${is_tarball}" == "true" ]; then
        # Delete the untarred database
        rm -rf database
    fi
    if [ "${task.ext.gtdb_keep_msa}" == "false" ]; then
        # Delete MSA of submitted and reference genomes.
        rm -rf supplemental/align/*.msa.fasta.gz
    fi
    rm -rf fna/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb-tk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
