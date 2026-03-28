#!/usr/bin/env nextflow
/**
 * Comprehensive screening of Klebsiella genomes for virulence and resistance determinants.
 *
 * This Bactopia Tool uses [Kleborate](https://github.com/katholt/Kleborate) to screen genome assemblies of
 * _Klebsiella pneumoniae_ and the _Klebsiella pneumoniae_ species complex (KpSC). Kleborate predicts:
 * MLST, species, ICEKp associated virulence loci, virulence plasmid associated loci,
 * antimicrobial resistance determinants, and K (capsule) and O antigen (LPS) serotype.
 *
 * @status stable
 * @keywords klebsiella, mlst, virulence, amr, serotyping, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,virulence,resistance,typing
 * @citation kleborate
 *
 * @subworkflows bactopiatool_init, kleborate
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Comprehensive Analysis
 * @publish *.kleborate.tsv      Comprehensive Kleborate report
 * @publish *-resistance.tsv     AMR determinant summary
 * @publish *-virulence.tsv      Virulence gene summary
 *
 * @section Merged Results
 * @publish kleborate.tsv        Merged TSV file containing Kleborate results from all samples
 *
 * @section Execution Logs
 * @publish logs/**               Tool execution logs
 * @publish logs/nf-*             Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml          Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { KLEBORATE         } from '../../../subworkflows/kleborate/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    KLEBORATE(BACTOPIATOOL_INIT.out.assembly)

    // Extract nf_logs as individual (meta, file) tuples for renaming
    ch_sample_nf_logs = KLEBORATE.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = KLEBORATE.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = KLEBORATE.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = KLEBORATE.out.run_outputs
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
