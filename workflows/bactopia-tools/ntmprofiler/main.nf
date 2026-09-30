#!/usr/bin/env nextflow
/**
 * Identify NTM species and predict antimicrobial resistance from sequencing reads.
 *
 * This Bactopia Tool uses [NTM-Profiler](https://github.com/jodyphelan/NTM-Profiler)
 * to identify nontuberculous Mycobacterium species and detect resistance-associated
 * variants using the default NTM-Profiler database or an optional custom database.
 *
 * @status stable
 * @keywords mycobacterium, ntm, species identification, antimicrobial resistance, surveillance
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation,database-dependent
 * @citation ntmprofiler
 *
 * @subworkflows utils_bactopia-tools, ntmprofiler
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input ntmprofiler_db?
 * Optional directory or compressed tarball created by `ntm-profiler update_db`.
 * When omitted, Bactopia downloads, prepares, and caches the database once.
 *
 * @section Per-Sample Results
 * @publish *.csv                      Per-sample profiling results in CSV format
 * @publish *.results.json.gz          Compressed machine-readable profiling results
 * @publish *.txt                      Human-readable profiling report
 * @publish supplemental/*             Alignments, variant calls, and other supporting files
 *
 * @section Merged Results
 * @publish ntmprofiler.csv             Collated species, lineage, and resistance results
 * @publish ntmprofiler.csv.variants.csv Collated variant calls
 * @publish ntmprofiler.csv.snp-dist.json Optional SNP-distance graph
 *
 * @section Execution Logs
 * @publish logs/ntmprofiler/*          Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                   Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                Software version information
 */
nextflow.enable.types = true

params {
    rundir : String
    ntmprofiler_db : Path?
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { NTMPROFILER         } from '../../../subworkflows/ntmprofiler/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_ntmprofiler = NTMPROFILER(ch_bactopiatool.reads, params.ntmprofiler_db)

    publish:
    // Per-sample
    sample_outputs = ch_ntmprofiler.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_ntmprofiler.sample_outputs)
    // Run-level
    run_outputs = ch_ntmprofiler.run_outputs
    run_nf_logs = collectNextflowLogs(ch_ntmprofiler.run_outputs)
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
