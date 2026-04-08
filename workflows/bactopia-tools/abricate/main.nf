#!/usr/bin/env nextflow
/**
 * Mass screening of contigs for antimicrobial resistance and virulence genes.
 *
 * This Bactopia Tool uses [Abricate](https://github.com/tseemann/abricate) to screen
 * assemblies against multiple resistance and virulence gene databases, including
 * NCBI, CARD, RESFINDER, ARG-ANNOT, VFDB, PLASMIDFINDER, ECOLI_VF, and MEGARES.
 * It processes a Bactopia analysis directory, runs Abricate on each sample, and
 * creates a merged summary report.
 *
 * @status stable
 * @keywords bacteria, antimicrobial resistance, virulence, screening, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation abricate, arg_annot, card, csvtk, ecoh, megares2, ncbi_reference_gene_catalog, plasmidfinder, resfinder, vfdb
 *
 * @subworkflows bactopiatool_init, abricate
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv            Tab-delimited report of Abricate screening results for each sample
 *
 * @section Merged Results
 * @publish abricate.tsv     Merged TSV report containing Abricate results from all samples
 *
 * @section Execution Logs
 * @publish logs/abricate/*  Tool execution logs (stdout/stderr)
 * @publish logs/nf-*        Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml     Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { ABRICATE            } from '../../../subworkflows/abricate/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_abricate = ABRICATE(ch_bactopiatool.assembly)

    publish:
    // Per-sample
    sample_outputs = ch_abricate.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_abricate.sample_outputs)
    // Run-level
    run_outputs = ch_abricate.run_outputs
    run_nf_logs = collectNextflowLogs(ch_abricate.run_outputs)
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
