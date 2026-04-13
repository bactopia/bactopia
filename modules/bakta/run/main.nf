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
 * @citation bakta, aragorn, diamond, hmmer, infernal, prodigal
 *
 * @note Database Required
 * Requires a Bakta database (directory or tarball) to be available.
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input db
 * Path to the Bakta database (Directory or compressed tarball)
 *
 * @input proteins?
 * FASTA file of trusted proteins to use for first-pass annotation
 *
 * @input prodigal_tf?
 * Prodigal training file for CDS prediction
 *
 * @input replicons?
 * Table (TSV/CSV) of replicon information for origin detection
 *
 * @output record(meta, blastdb, faa, ffn, fna, gbff, gff, hypotheticals_tsv, hypotheticals_faa, inference_tsv, json, png, svg, tsv, txt, results, logs, nf_logs, versions)
 * - `blastdb`: A compressed tar.gz archive of BLAST+ databases of the contigs, genes, and proteins
 * - `faa`: CDS/sORF amino acid sequences as FASTA
 * - `ffn`: Feature nucleotide sequences as FASTA
 * - `fna`: Replicon/contig DNA sequences as FASTA
 * - `gbff`: Annotations and sequences in GenBank format
 * - `gff`: Annotations and sequences in GFF3 format
 * - `hypotheticals_tsv`: Further information on hypothetical protein CDS as tab-separated values
 * - `hypotheticals_faa`: Hypothetical protein CDS amino acid sequences as FASTA
 * - `inference_tsv`: Detailed annotation evidence and database hit information
 * - `json`: Machine-readable annotations and metadata in JSON format
 * - `png`: Circular genome plot as PNG image
 * - `svg`: Circular genome plot as SVG image
 * - `tsv`: Annotations as simple human readable tab-separated values
 * - `txt`: Broad summary of Bakta annotations
 *
 * @results additional
 * - `${prefix}.embl`: Annotations and sequences in EMBL format
 * - `${prefix}.inference.tsv`: Detailed annotation evidence and database hit information
 * - `${prefix}.json`: Machine-readable annotations and metadata in JSON format
 * - `${prefix}.png`: Circular genome plot as PNG image
 * - `${prefix}.svg`: Circular genome plot as SVG image
 */
nextflow.preview.types = true

process BAKTA_RUN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )
    db         : Path
    proteins   : Path?
    prodigal_tf: Path?
    replicons  : Path?

    stage:
    stageAs fna, 'staging/fna/*'

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        blastdb: file("${prefix}-blastdb.tar.gz"),
        faa: file("${prefix}.{faa,faa.gz}"),
        ffn: file("${prefix}.{ffn,ffn.gz}"),
        fna: file("${prefix}.{fna,fna.gz}"),
        gbff: file("${prefix}.{gbff,gbff.gz}"),
        gff: file("${prefix}.{gff3,gff3.gz}"),
        hypotheticals_tsv: file("${prefix}.hypotheticals.tsv"),
        hypotheticals_faa: file("${prefix}.hypotheticals.{faa,faa.gz}"),
        inference_tsv: file("${prefix}.inference.tsv"),
        json: file("${prefix}.{json,json.gz}"),
        png: file("${prefix}.png"),
        svg: file("${prefix}.{svg,svg.gz}"),
        tsv: file("${prefix}.tsv"),
        txt: file("${prefix}.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}-blastdb.tar.gz"),
            files("${prefix}.{embl,embl.gz}"),
            files("${prefix}.{faa,faa.gz}"),
            files("${prefix}.{ffn,ffn.gz}"),
            files("${prefix}.{fna,fna.gz}"),
            files("${prefix}.{gbff,gbff.gz}"),
            files("${prefix}.{gff3,gff3.gz}"),
            files("${prefix}.hypotheticals.tsv"),
            files("${prefix}.hypotheticals.{faa,faa.gz}"),
            files("${prefix}.inference.tsv"),
            files("${prefix}.{json,json.gz}"),
            files("${prefix}.png"),
            files("${prefix}.{svg,svg.gz}"),
            files("${prefix}.tsv"),
            files("${prefix}.txt")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/main/annotator/bakta/",
        logs_dir: "${prefix}/main/annotator/bakta/logs/",
        process_name: task.ext.process_name
    )

    def proteins_opt = proteins != null ? "--proteins ${proteins.getName()}" : ""
    def prodigal_opt = prodigal_tf != null ? "--prodigal-tf ${prodigal_tf.getName()}" : ""
    def replicons_opt = replicons != null ? "--replicons ${replicons.getName()}" : ""
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
        ${fna}

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
        gzip --best bakta/${prefix}.json
        gzip --best bakta/${prefix}.svg
    fi

    # Cleanup
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database/
    fi
    rm -rf blastdb/
    mv bakta/*.log ./
    mv bakta/* ./
    rmdir bakta/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
