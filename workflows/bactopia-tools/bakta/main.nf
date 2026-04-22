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
 * @subworkflows utils_bactopia-tools, bakta
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
 * @section Execution Logs
 * @publish logs/**                Tool execution logs
 * @publish logs/nf-*              Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
   */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    bakta_db              : Value<Path>
    download_bakta        : Boolean
    bakta_save_as_tarball : Boolean
    bakta_proteins        : Value<Path?>
    bakta_prodigal_tf     : Value<Path?>
    replicons             : Value<Path?>
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { BAKTA               } from '../../../subworkflows/bakta/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_bakta = BAKTA(
        ch_bactopiatool.assembly,
        params.bakta_db,
        params.download_bakta,
        params.bakta_save_as_tarball,
        params.bakta_proteins,
        params.bakta_prodigal_tf,
        params.replicons
    )

    publish:
    // Per-sample
    sample_outputs = ch_bakta.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_bakta.sample_outputs)
    // Run-level
    run_outputs = ch_bakta.run_outputs
    run_nf_logs = collectNextflowLogs(ch_bakta.run_outputs)
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
