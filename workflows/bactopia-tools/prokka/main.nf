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
 * @subworkflows utils_bactopia-tools, prokka
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input prokka_proteins
 * Path to trusted protein FASTA file for additional homology-based annotation
 *
 * @input prokka_prodigal_tf
 * Path to a Prodigal training file for gene prediction
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
    rundir : String

    // Tool-specific parameters
    prokka_proteins    : Path?
    prokka_prodigal_tf : Path?
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { PROKKA              } from '../../../subworkflows/prokka/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_prokka = PROKKA(
        ch_bactopiatool.assembly,
        params.prokka_proteins,
        params.prokka_prodigal_tf
    )

    publish:
    // Per-sample
    sample_outputs = ch_prokka.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_prokka.sample_outputs)
    // Run-level
    run_outputs = ch_prokka.run_outputs
    run_nf_logs = collectNextflowLogs(ch_prokka.run_outputs)
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results.flatten()  >> "${r.meta.output_dir}/"
            r.logs.flatten()     >> "${r.meta.logs_dir}/"
            r.versions.flatten() >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results.flatten()  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
