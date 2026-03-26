#!/usr/bin/env nextflow
/**
 * Rapid annotation of bacterial genomes and plasmids.
 *
 * This Bactopia Tool uses [Bakta](https://github.com/oschwengers/bakta) to rapidly annotate bacterial
 * genomes and plasmids in a standardized fashion. Bakta makes use of a large database ([40+ GB](https://doi.org/10.5281/zenodo.4247252))
 * to provide extensive annotations including: tRNA, tmRNA, rRNA, ncRNA, CRISPR, CDS, and sORFs.
 *
 * @status stable
 * @keywords bacteria, fasta, annotation, genbank, gff, proteins, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,annotation,database-dependent
 * @citation bakta
 *
 * @subworkflows bactopiatool_init, bakta
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input bakta_db
 * Path to Bakta database for genome annotation
 *
 * @input download_bakta
 * Download Bakta database if not found locally
 *
 * @input bakta_save_as_tarball
 * Save Bakta database as compressed tarball for reuse
 *
 * @input bakta_proteins
 * Additional protein sequences for homology search
 *
 * @input bakta_prodigal_tf
 * Path to Prodigal training file for gene prediction
 *
 * @input replicons
 * Additional replicon sequences for contamination screening
 *
 * @section Annotation
 * @publish *.gff3                 Genome annotation in GFF3 format
 * @publish *.gbff                 Genome annotation in GenBank format
 * @publish *.faa                  Protein sequences in FASTA format
 * @publish *.ffn                  Feature nucleotide sequences
 * @publish *.fna                  Nucleotide sequences of all features
 * @publish *.hypotheticals.tsv    List of hypothetical proteins
 * @publish *.tsv                  Annotation summary in TSV format
 * @publish *.txt                  Detailed annotation report
 *
 * @section Merged Results
 * @note Merged results are not created for Bakta annotations
 *
 * @section Execution Logs
 * @publish logs/**                Tool execution logs
 * @publish logs/nf-*              Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
   */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    bakta_db              : Path
    download_bakta        : Boolean
    bakta_save_as_tarball : Boolean
    bakta_proteins        : Path?
    bakta_prodigal_tf     : Path?
    replicons             : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { BAKTA             } from '../../../subworkflows/bakta/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    BAKTA(
        BACTOPIATOOL_INIT.out.assembly,
        params.bakta_db,
        params.download_bakta,
        params.bakta_save_as_tarball,
        params.bakta_proteins,
        params.bakta_prodigal_tf,
        params.replicons
    )

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = BAKTA.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = BAKTA.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = BAKTA.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = BAKTA.out.run_outputs
    run_nf_logs = ch_run_nf_logs
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
