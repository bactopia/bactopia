#!/usr/bin/env nextflow
/**
 * Detection of antimicrobial resistance and lineage typing of Mycobacterium tuberculosis.
 *
 * This Bactopia Tool uses [TBProfiler](https://github.com/jodyphelan/TBProfiler) to profile
 * Mycobacterium tuberculosis genomes for resistance mutations and strain typing. The workflow
 * processes sequencing reads to identify resistance-conferring variants and determine
 * the lineage of each TB isolate.
 *
 * @status stable
 * @keywords Mycobacterium tuberculosis, resistance, lineage, typing, TB
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation tbprofiler
 *
 * @subworkflows bactopiatool_init, tbprofiler
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.results.txt              Text file containing TBProfiler resistance and strain typing results
 * @publish *.results.json            JSON file containing detailed TBProfiler analysis results
 * @publish *.results.csv             CSV file containing TBProfiler results in tabular format
 * @publish bam/*.bam                 BAM file with read alignment details against reference genomes
 * @publish vcf/*.targets.csq.vcf.gz   VCF file with variant annotations and functional consequences
 *
 * @section Merged Results
 * @publish tbprofiler.tsv             Merged TSV file containing TBProfiler results from all samples
 *
 * @section Execution Logs
 * @publish logs/tbprofiler/*          Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { TBPROFILER        } from '../../../subworkflows/tbprofiler/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    TBPROFILER(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(TBPROFILER.out.results)
    ch_logs = ch_logs.mix(TBPROFILER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(TBPROFILER.out.nf_logs)
    ch_versions = ch_versions.mix(TBPROFILER.out.versions)

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
