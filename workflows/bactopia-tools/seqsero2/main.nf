#!/usr/bin/env nextflow
/**
 * Salmonella serotype prediction from sequencing reads or assemblies.
 *
 * This Bactopia Tool uses [SeqSero2](https://github.com/denglab/SeqSero2) to predict Salmonella
 * serotypes from both raw sequencing reads and assembled genomes. SeqSero2 is a novel pipeline
 * for determining Salmonella serotypes using raw sequencing reads or assemblies through
 * k-mer analysis and targeted identification of O and H antigen genes. The tool provides
 * accurate serotype predictions following the Kaufmann-White scheme, supporting traditional
 * and molecular serotyping methods for epidemiological surveillance and outbreak investigation.
 *
 * @status stable
 * @keywords Salmonella, serotyping, epidemiology, O antigen, H antigen, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, seqsero2
 *
 * @subworkflows bactopiatool_init, seqsero2
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *_result.tsv             Tab-delimited file with detailed SeqSero2 results for each sample
 * @publish *_result.txt             Text file with key-value pairs of SeqSero2 prediction results
 * @publish *_log.txt                Detailed log file from SeqSero2 analysis
 *
 * @section Merged Results
 * @publish seqsero2.tsv             Merged TSV file containing SeqSero2 results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                  Tool execution logs including SeqSero2 output
 * @publish logs/nf-*                Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml             Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SEQSERO2          } from '../../../subworkflows/seqsero2/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    SEQSERO2(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(SEQSERO2.out.results)
    ch_logs = ch_logs.mix(SEQSERO2.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SEQSERO2.out.nf_logs)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions)

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
