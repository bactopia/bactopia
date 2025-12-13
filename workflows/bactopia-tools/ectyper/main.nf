#!/usr/bin/env nextflow
/**
 * In silico prediction of Escherichia coli serotype.
 *
 * This Bactopia Tool uses [ECTyper](https://github.com/phac-nml/ecoli_serotyping) to conduct
 * _in silico_ prediction of serotype for _Escherichia coli_ genomes. It uses the genome assemblies
 * to provide basic species identification and the predicted _E. coli_ serotype (e.g. O174:H21).
 *
 * @status stable
 * @keywords escherichia coli, serotyping, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation ectyper
 *
 * @subworkflows bactopiatool_init, ectyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                  Tab-delimited file with ECTyper result
 * @publish blast_output_alleles.txt Allele report generated from BLAST results
 *
 * @section Merged Results
 * @publish ectyper.tsv            Merged TSV file containing ECTyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                 Tool execution logs
 * @publish logs/nf-*               Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml            Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { ECTYPER           } from '../../../subworkflows/ectyper/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    ECTYPER(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(ECTYPER.out.results)
    ch_logs = ch_logs.mix(ECTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(ECTYPER.out.nf_logs)
    ch_versions = ch_versions.mix(ECTYPER.out.versions)

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
