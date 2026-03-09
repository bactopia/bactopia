#!/usr/bin/env nextflow
/**
 * A NATA accredited tool for reporting the presence of antimicrobial resistance genes.
 *
 * This Bactopia Tool uses [abriTAMR](https://github.com/MDU-PHL/abritamr) to identify
 * antimicrobial resistance genes in bacterial genomes. It runs AMRFinderPlus on each
 * sample and collates the results into functional classes, producing detailed reports
 * on resistance genes, partial matches, and virulence factors. It is accredited by NATA
 * for use in reporting presence of reportable AMR genes in Victoria, Australia.
 *
 * @status stable
 * @keywords bacteria, antimicrobial resistance, virulence, amr, nata, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation abritamr, amrfinderplus, csvtk
 *
 * @subworkflows bactopiatool_init, abritamr
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.abritamr.txt           Tab-delimited file combining non-empty summary files from abriTAMR
 * @publish *.amrfinder.out          Raw output from AMRFinderPlus (per sequence)
 * @publish *.summary_matches.txt    Tab-delimited file with AMR gene matches per sequence
 * @publish *.summary_partials.txt   Tab-delimited file with partial AMR gene matches
 * @publish *.summary_virulence.txt  Tab-delimited file with virulence gene classifications
 *
 * @section Merged Results
 * @publish abritamr.tsv             Merged TSV file containing AMR summaries from all samples
 *
 * @section Execution Logs
 * @publish logs/abritamr/*          Tool execution logs (stdout/stderr)
 * @publish logs/nf-*                Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml             Software version information
 */
nextflow.preview.types = true

params {
    rundir: String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { ABRITAMR          } from '../../../subworkflows/abritamr/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    ABRITAMR(BACTOPIATOOL_INIT.out.assembly)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = ABRITAMR.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = ABRITAMR.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ABRITAMR.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = ABRITAMR.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.results  >> "${r.meta.output_dir}/"
            r.logs     >> "${r.meta.logs_dir}/"
            r.versions >> "${r.meta.logs_dir}/"
        }
    }
    sample_nf_logs {
        path { meta, f -> f >> "${meta.logs_dir}/nf${f.name}" }
    }

    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_outputs {
        path { r ->
            r.results  >> "${params.rundir}/${r.meta.output_dir}/"
            r.logs     >> "${params.rundir}/${r.meta.logs_dir}/"
            r.versions >> "${params.rundir}/${r.meta.logs_dir}/"
        }
    }
    run_nf_logs {
        path { meta, f -> f >> "${params.rundir}/${meta.logs_dir}/nf${f.name}" }
    }
}
