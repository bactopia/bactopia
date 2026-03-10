/**
 * Rapid and standardized annotation of bacterial genomes and plasmids.
 *
 * Uses [Bakta](https://github.com/oschwengers/bakta) to annotate genomes via alignment-free
 * sequence identification. It detects CDS, sORFs, tRNAs, tmRNAs, rRNAs, ncRNAs, and CRISPR
 * arrays, assigning functions from a comprehensive database.
 *
 * @status stable
 * @keywords bacteria, annotation, genome, assembly, prodigal, compliant, genbank, ena
 * @tags complexity:complex input-type:multiple output-type:multiple features:database-dependent,conditional-logic,archive-output
 * @citation bakta
 *
 * @note Database Required
 * Requires a Bakta database (directory or tarball) to be available.
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * Path to the Bakta database (Directory or compressed tarball)
 *
 * @input proteins
 * Optional FASTA file of trusted proteins to use for first-pass annotation
 *
 * @input prodigal_tf
 * Optional Prodigal training file for CDS prediction
 *
 * @input replicons
 * Optional table (TSV/CSV) of replicon information for origin detection
 *
 * @output record(meta, annotations, embl, faa, ffn, fna, gbff, gff, hypotheticals_tsv, hypotheticals_faa, tsv, txt, blastdb, results, logs, nf_logs, versions)
 * - `embl`: Annotations and sequences in EMBL format
 * - `faa`: CDS/sORF amino acid sequences as FASTA
 * - `ffn`: Feature nucleotide sequences as FASTA
 * - `fna`: Replicon/contig DNA sequences as FASTA
 * - `gbff`: Annotations and sequences in GenBank format
 * - `gff`: Annotations and sequences in GFF3 format
 * - `hypotheticals_tsv`: Further information on hypothetical protein CDS as tab-separated values
 * - `hypotheticals_faa`: Hypothetical protein CDS amino acid sequences as FASTA
 * - `tsv`: Annotations as simple human readable tab-separated values
 * - `txt`: Broad summary of Bakta annotations
 * - `blastdb`: A compressed tar.gz archive of BLAST+ databases of the contigs, genes, and proteins
 */
nextflow.preview.types = true

process BAKTA_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, assembly: Path): Record
    db             : Path
    proteins       : Path?
    prodigal_tf    : Path?
    replicons      : Path?

    output:
    record(
        meta: meta,
        annotations: [file("bakta/${prefix}.{fna,fna.gz}"), file("bakta/${prefix}.{faa,faa.gz}"), file("bakta/${prefix}.{gff3,gff3.gz}")],
        embl: file("bakta/${prefix}.{embl,embl.gz}"),
        faa: file("bakta/${prefix}.{faa,faa.gz}"),
        ffn: file("bakta/${prefix}.{ffn,ffn.gz}"),
        fna: file("bakta/${prefix}.{fna,fna.gz}"),
        gbff: file("bakta/${prefix}.{gbff,gbff.gz}"),
        gff: file("bakta/${prefix}.{gff3,gff3.gz}"),
        hypotheticals_tsv: file("bakta/${prefix}.hypotheticals.tsv"),
        hypotheticals_faa: file("bakta/${prefix}.hypotheticals.{faa,faa.gz}"),
        tsv: file("bakta/${prefix}.tsv"),
        txt: file("bakta/${prefix}.txt"),
        blastdb: file("bakta/${prefix}-blastdb.tar.gz"),
        results: [file("bakta/${prefix}.tsv"), file("bakta/${prefix}.txt"), file("bakta/${prefix}.{embl,embl.gz}"), file("bakta/${prefix}.{faa,faa.gz}"), file("bakta/${prefix}.{ffn,ffn.gz}"), file("bakta/${prefix}.{fna,fna.gz}"), file("bakta/${prefix}.{gbff,gbff.gz}"), file("bakta/${prefix}.{gff3,gff3.gz}"), file("bakta/${prefix}.hypotheticals.tsv"), file("bakta/${prefix}.hypotheticals.{faa,faa.gz}"), file("bakta/${prefix}-blastdb.tar.gz")],
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
    meta.output_dir = "${prefix}/main/annotator/bakta/"
    meta.logs_dir = "${prefix}/main/annotator/bakta/logs/"
    meta.process_name = task.ext.process_name

    def proteins_opt = proteins.getName() != "EMPTY_PROTEINS" ? "--proteins ${proteins.getName()}" : ""             // TODO: Remove when Path? is fixed
    def prodigal_opt = prodigal_tf.getName() != "EMPTY_PRODIGAL_TF" ? "--prodigal-tf ${prodigal_tf.getName()}" : "" // TODO: Remove when Path? is fixed
    def replicons_opt = replicons.getName() != "EMPTY_REPLICONS" ? "--replicons ${replicons.getName()}" : ""        // TODO: Remove when Path? is fixed
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
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
        ${assembly}

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
