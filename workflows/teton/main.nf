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
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIA_INIT()

    // Gather samples in one place
    GATHER(BACTOPIA_INIT.out.samples)
    ch_results = ch_results.mix(GATHER.out.results)
    ch_logs = ch_logs.mix(GATHER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GATHER.out.nf_logs)
    ch_versions = ch_versions.mix(GATHER.out.versions)

    // Run Teton
    TETON(
        GATHER.out.fastq_only,
        params.kraken2_db,
        params.use_srascrubber
    )
    ch_results = ch_results.mix(TETON.out.results)
    ch_logs = ch_logs.mix(TETON.out.logs)
    ch_nf_logs = ch_nf_logs.mix(TETON.out.nf_logs)
    ch_versions = ch_versions.mix(TETON.out.versions)

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
