#!/usr/bin/env nextflow
/**
 * Assessment of microbial genome assembly quality.
 *
 * This Bactopia Tool uses [CheckM](https://github.com/Ecogenomics/CheckM) to assess the quality
 * of microbial genomes recovered from isolates, single cells, and metagenomes using
 * a set of lineage-specific marker genes.
 *
 * @status stable
 * @keywords assembly quality, microbial genomes, completeness, contamination, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,assembly-qa
 * @citation checkm
 *
 * @subworkflows bactopiatool_init, checkm
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Quality Assessment
 * @publish *.genes.aln          Alignment of multi-copy genes and their AAI identity
 * @publish *.results.txt        Final results of CheckM's lineage_wf
 * @publish lineage.ms           Output file describing marker set for each bin
 * @publish bins/**              Directory with inputs for processing by CheckM
 * @publish storage/**           Directory with intermediate results from CheckM processing
 *
 * @section Merged Results
 * @publish checkm.tsv           Merged TSV file with CheckM results from all samples
 *
 * @section Execution Logs
 * @publish logs/**              Tool execution logs
 * @publish logs/nf-*            Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml         Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { CHECKM            } from '../../../subworkflows/checkm/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    CHECKM(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(CHECKM.out.results)
    ch_logs = ch_logs.mix(CHECKM.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CHECKM.out.nf_logs)
    ch_versions = ch_versions.mix(CHECKM.out.versions)

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
