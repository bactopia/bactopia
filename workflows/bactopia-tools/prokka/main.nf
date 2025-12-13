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
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    PROKKA(
        BACTOPIATOOL_INIT.out.samples,
        params.prokka_proteins,
        params.prokka_prodigal_tf
    )

    // Collect outputs
    ch_results = ch_results.mix(PROKKA.out.results)
    ch_logs = ch_logs.mix(PROKKA.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PROKKA.out.nf_logs)
    ch_versions = ch_versions.mix(PROKKA.out.versions)

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}
