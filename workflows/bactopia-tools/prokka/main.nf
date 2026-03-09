#!/usr/bin/env nextflow
/**
 * Rapid whole genome annotation of bacterial, archaeal, and viral genomes.
 *
 * This Bactopia Tool uses [Prokka](https://github.com/tseemann/prokka) to rapidly annotate small genomes
 * in a standardized fashion. It identifies protein-coding genes, rRNA, tRNA, and other features,
 * then searches them against multiple reference databases to provide comprehensive functional annotation.
 *
 * @status stable
 * @keywords annotation, genome, prokaryote, functional annotation, genes
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool
 * @citation prokka
 *
 * @subworkflows bactopiatool_init, prokka
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.gff                 Genome annotation in GFF3 format containing sequences and annotations
 * @publish *.gbk                 Genome annotation in GenBank format containing sequences and annotations
 * @publish *.faa                 Protein FASTA file of translated CDS sequences
 * @publish *.fna                 Nucleotide FASTA file of input contig sequences
 * @publish *.ffn                 Nucleotide FASTA file of all predicted transcripts
 * @publish *.fsa                 Nucleotide FASTA file of predicted protein sequences
 * @publish *.sqn                 ASN1 format Sequin file for GenBank submission
 * @publish *.tbl                 Feature Table file for GenBank submission
 * @publish *.tsv                 Tab-separated file of all features with functional information
 * @publish *.txt                 Statistics report of annotated features
 * @publish *.blastdb.tar.gz      BLAST+ database archive of contigs, genes, and proteins
 *
 * @section Merged Results
 * @publish prokka.tsv            Merged TSV file containing annotation summaries from all samples
 *
 * @section Execution Logs
 * @publish logs/prokka/*         Tool execution logs (stdout/stderr)
 * @publish logs/nf-*             Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml          Software version information
 */
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    prokka_proteins    : Path?
    prokka_prodigal_tf : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { PROKKA            } from '../../../subworkflows/prokka/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    PROKKA(
        BACTOPIATOOL_INIT.out.assembly,
        params.prokka_proteins,
        params.prokka_prodigal_tf
    )
    ch_sample_nf_logs = PROKKA.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    sample_outputs = PROKKA.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
}

output {
    sample_outputs {
        path { r ->
            r.results  >> "${r.meta.output_dir}/"
            r.logs     >> "${r.meta.logs_dir}/"
            r.versions >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }
}
