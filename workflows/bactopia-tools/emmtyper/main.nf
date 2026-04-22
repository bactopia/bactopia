#!/usr/bin/env nextflow
/**
 * emm-typing of Streptococcus pyogenes assemblies.
 *
 * This Bactopia Tool uses [emmtyper](https://github.com/MDU-PHL/emmtyper) for
 * emm-typing of _Streptococcus pyogenes_ using a de novo or complete assembly.
 *
 * @status stable
 * @keywords streptococcus pyogenes, emm typing, fasta, bactopia-tool
 * @tags complexity:simple input-type:parameter output-type:multiple features:bactopia-tool,typing
 * @citation emmtyper
 *
 * @subworkflows utils_bactopia-tools, emmtyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input emmtyper_blastdb
 * Path to emmtyper BLAST database
 *
 * @section Per-Sample Results
 * @publish *.tsv           Tab-delimited file with emmtyper result
 *
 * @section Merged Results
 * @publish emmtyper.tsv     Merged TSV file containing emmtyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**          Tool execution logs
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    emmtyper_blastdb : Value<Path?>
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { EMMTYPER            } from '../../../subworkflows/emmtyper/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_emmtyper = EMMTYPER(ch_bactopiatool.assembly, params.emmtyper_blastdb)

    publish:
    // Per-sample
    sample_outputs = ch_emmtyper.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_emmtyper.sample_outputs)
    // Run-level
    run_outputs = ch_emmtyper.run_outputs
    run_nf_logs = collectNextflowLogs(ch_emmtyper.run_outputs)
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
