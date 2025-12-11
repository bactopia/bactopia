/**
 * Bakta is a rapid & standardized annotation of bacterial genomes & plasmids.
 *
 * This process executes bakta_run to perform analysis
 *
 * @status stable
 * @keywords annotation, fasta, genome
 * @tags complexity:complex input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, database-dependent, path-workarounds
 * @citation bakta_run
 *
 * @note Uses EMPTY_* placeholder files for optional parameters
 * @note Requires external database to be available
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: Genome assembly in FASTA format
 *
 * @input db
 * Bakta database or tarball
 *
 * @input proteins
 * Fasta file of trusted proteins to first annotate from
 *
 * @input prodigal_tf
 * Training file to use for CDS prediction
 *
 * @input replicons
 * Replicons/plasmids file for replicon detection
 *
 * @output annotations       Tuple containing the main annotation files (fna, faa, gff3)
 * @output embl              Annotations in EMBL format
 * @output faa               Protein CDS sequences
 * @output ffn               Nucleotide CDS sequences
 * @output fna               Nucleotide sequences of the features
 * @output gbff              Annotations in GenBank format
 * @output gff               Annotations in GFF3 format
 * @output hypotheticals_tsv Summary of hypothetical proteins
 * @output hypotheticals_faa Hypothetical protein sequences
 * @output tsv               Summary of annotated features
 * @output txt               Sequence and annotation statistics
 * @output blastdb           BLAST database of contigs, genes, and proteins
 * @output logs              Optional tool execution logs
 * @output nf_logs           Nextflow execution logs
 * @output versions          Software version information (YAML format)
 */
nextflow.preview.types = true

process BAKTA_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Set<Path>>
    db             : Set<Path>?
    proteins       : Path?
    prodigal_tf    : Path?
    replicons      : Path?

    output:
    annotations       = tuple(meta, files("bakta/${prefix}.{fna,fna.gz}"), files("bakta/${prefix}.{faa,faa.gz}"), files("bakta/${prefix}.{gff3,gff3.gz}"))
    embl              = tuple(meta, file("bakta/${prefix}.{embl,embl.gz}"))
    faa               = tuple(meta, file("bakta/${prefix}.{faa,faa.gz}"))
    ffn               = tuple(meta, file("bakta/${prefix}.{ffn,ffn.gz}"))
    fna               = tuple(meta, file("bakta/${prefix}.{fna,fna.gz}"))
    gbff              = tuple(meta, file("bakta/${prefix}.{gbff,gbff.gz}"))
    gff               = tuple(meta, file("bakta/${prefix}.{gff3,gff3.gz}"))
    hypotheticals_tsv = tuple(meta, file("bakta/${prefix}.hypotheticals.tsv"))
    hypotheticals_faa = tuple(meta, file("bakta/${prefix}.hypotheticals.{faa,faa.gz}"))
    tsv               = tuple(meta, file("bakta/${prefix}.tsv"))
    txt               = tuple(meta, file("bakta/${prefix}.txt"))
    blastdb           = tuple(meta, files("bakta/${prefix}-blastdb.tar.gz"))
    logs              = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs           = tuple(meta, files(".command.*"))
    versions          = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/main/annotator/bakta/"
    meta.logs_dir = "${prefix}/main/annotator/bakta/logs/"
    meta.process_name = task.ext.process_name
    def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""             // TODO: Remove when Path? is fixed
    def prodigal_opt = prodigal_tf.getName() != "EMPTY_PRODIGAL_TF" ? "--prodigal-tf ${prodigal_tf.getName()}" : "" // TODO: Remove when Path? is fixed
    def replicons_opt = replicons.getName() != "EMPTY_REPLICONS" ? "--replicons ${replicons.getName()}" : ""        // TODO: Remove when Path? is fixed
    def is_tarball = db.toList()[0].getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        BAKTA_DB=\$(find database/ -name "bakta.db" | sed 's=bakta\\.db==')
    else
        BAKTA_DB=\$(find ${db}/ -name "bakta.db" | sed 's=bakta\\.db==')
    fi

    bakta \\
        --output bakta \\
        ${task.ext.args} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        --db \$BAKTA_DB \\
        ${proteins_opt} \\
        ${prodigal_opt} \\
        ${replicons_opt} \\
        ${fasta}

    # Make blastdb of contigs, genes, proteins
    mkdir blastdb
    cat bakta/${prefix}.fna | makeblastdb -dbtype "nucl" -title "Assembled contigs for ${prefix}" -out blastdb/${prefix}.fna
    cat bakta/${prefix}.ffn | makeblastdb -dbtype "nucl" -title "Predicted genes sequences for ${prefix}" -out blastdb/${prefix}.ffn
    cat bakta/${prefix}.faa | makeblastdb -dbtype "prot" -title "Predicted protein sequences for ${prefix}" -out blastdb/${prefix}.faa
    tar -cvf - blastdb/ | gzip -c > bakta/${prefix}-blastdb.tar.gz

    if [[ "${task.ext.skip_compression}" == "false" ]]; then
        gzip --best bakta/${prefix}.embl
        gzip --best bakta/${prefix}.faa
        gzip --best bakta/${prefix}.ffn
        gzip --best bakta/${prefix}.fna
        gzip --best bakta/${prefix}.gbff
        gzip --best bakta/${prefix}.gff3
        gzip --best bakta/${prefix}.hypotheticals.faa
    fi

    # Clean up
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database
    fi
    rm -rf blastdb/
    mv bakta/*.log ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
