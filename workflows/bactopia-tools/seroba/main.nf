#!/usr/bin/env nextflow
/**
 * Serotyping of Streptococcus pneumoniae from Illumina paired-end reads.
 *
 * This Bactopia Tool uses [Seroba](https://github.com/sanger-pathogens/seroba) to predict the
 * serotype of *Streptococcus pneumoniae* samples from raw sequencing reads. Seroba employs
 * a k-mer based approach to identify and type pneumococcal capsules, determining both the
 * serotype and serogroup based on the presence of specific capsular polysaccharide synthesis
 * (cps) locus sequences. The tool is specifically designed for Illumina paired-end reads and
 * provides accurate serotype predictions essential for epidemiological surveillance and vaccine
 * development studies.
 *
 * @status stable
 * @keywords serotyping, Streptococcus pneumoniae, capsule, cps locus, vaccine, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, seroba
 *
 * @subworkflows bactopiatool_init, seroba
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                              Tab-delimited file containing predicted serotype for each sample
 * @publish *detailed_serogroup_info.txt      Detailed information about serotype prediction and coverage
 *
 * @section Merged Results
 * @publish seroba.tsv                         Merged TSV file containing serotype predictions from all samples
 *
 * @section Execution Logs
 * @publish logs/**                            Tool execution logs including Seroba output
 * @publish logs/nf-*                          Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                       Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SEROBA            } from '../../../subworkflows/seroba/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    SEROBA(BACTOPIATOOL_INIT.out.assembly)

    // Collect outputs
    ch_results = ch_results.mix(SEROBA.out.results)
    ch_logs = ch_logs.mix(SEROBA.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SEROBA.out.nf_logs)
    ch_versions = ch_versions.mix(SEROBA.out.versions)

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
