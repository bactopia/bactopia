#!/usr/bin/env nextflow
/**
 * Taxonomic classification and abundance profiling of metagenomic reads.
 *
 * This workflow performs metagenomic classification using [Kraken2](https://github.com/DerrickWood/kraken2)
 * and [Bracken](https://github.com/jenniferlu717/Bracken), with optional host read removal
 * using SRA Scrubber. It processes metagenomic sequencing reads to estimate bacterial
 * genome sizes and separate bacterial from non-bacterial organisms.
 *
 * @status stable
 * @keywords metagenomics, classification, kraken2, bracken, abundance, profiling
 * @tags complexity:complex input-type:parameter output-type:multiple features:aggregation,conditional-logic,database-dependent
 *
 * @subworkflows bactopia_init, gather, teton
 *
 * @input rundir
 * Directory containing metagenomic sequencing reads
 *
 * @input kraken2_db
 * Path to Kraken2 database for classification
 *
 * @input use_srascrubber
 * Remove host reads using SRA scrubber before classification
 *
 * @section Per-Sample Results
 * @publish bacteria.tsv               Per-sample TSV files containing bacterial organisms and their properties
 * @publish nonbacteria.tsv            Per-sample TSV files containing non-bacterial organisms
 * @publish sizemeup.tsv               Per-sample TSV files with genome size estimates
 *
 * @section Merged Results
 * @publish merged-bacteria.tsv        Consolidated TSV file of all bacterial organisms across samples
 * @publish merged-nonbacteria.tsv     Consolidated TSV file of all non-bacterial organisms across samples
 * @publish merged-sizemeup.tsv        Consolidated TSV file of genome size estimates across samples
 * @publish report.tsv                 Joined TSV file combining scrubber and classification results
 *
 * @section Execution Logs
 * @publish logs/**                    Tool execution logs
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    rundir   : String

    // Tool-specific parameters
    kraken2_db      : Path
    use_srascrubber : Boolean
}

include { BACTOPIA_INIT   } from '../../subworkflows/utils/bactopia'
include { GATHER          } from '../../subworkflows/bactopia/gather/main'
include { TETON           } from '../../subworkflows/teton/main'

workflow {
    main:
    BACTOPIA_INIT()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)

    // Run Teton
    TETON(
        GATHER.out.reads,
        params.kraken2_db,
        params.use_srascrubber,
        params.nohuman_db ? file(params.nohuman_db) : file("NO_DB"),
        params.download_nohuman,
        params.nohuman_save_as_tarball
    )

    // Collect all sample-level outputs
    ch_sample_outputs = GATHER.out.sample_outputs
        .mix(TETON.out.scrubber_outputs)
        .mix(TETON.out.bracken_outputs)
        .mix(TETON.out.samplesheet_outputs)
        .mix(TETON.out.report_outputs)

    // Collect all run-level outputs
    ch_run_outputs = GATHER.out.run_outputs
        .mix(TETON.out.merged_report)
        .mix(TETON.out.merged_bacteria)
        .mix(TETON.out.merged_nonbacteria)
        .mix(TETON.out.merged_sizemeup)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = ch_sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = ch_run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ch_sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = ch_run_outputs
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
