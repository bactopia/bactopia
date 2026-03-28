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
    BACTOPIATOOL_INIT()
    CHECKM(BACTOPIATOOL_INIT.out.assembly)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = CHECKM.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = CHECKM.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = CHECKM.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = CHECKM.out.run_outputs
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
