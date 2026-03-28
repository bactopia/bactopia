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
    BACTOPIATOOL_INIT()
    TBPROFILER(BACTOPIATOOL_INIT.out.reads)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = TBPROFILER.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = TBPROFILER.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = TBPROFILER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = TBPROFILER.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.bam      >> "${r.meta.output_dir}/"
            r.csv      >> "${r.meta.output_dir}/"
            r.json     >> "${r.meta.output_dir}/"
            r.txt      >> "${r.meta.output_dir}/"
            r.vcf      >> "${r.meta.output_dir}/"
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
            r.itol     >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs.flatten()     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions.flatten() >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
