#!/usr/bin/env nextflow
/**
 * Rapid determination of Shigella serotypes from sequencing reads.
 *
 * This Bactopia Tool uses [ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper) to rapidly
 * determine *Shigella* serotypes from both Illumina (single or paired-end) and Oxford Nanopore
 * reads. ShigaTyper performs k-mer based analysis targeting specific antigenic determinants
 * and marker genes to predict serotypes according to the White-Kauffmann-Le Minor classification
 * scheme. The tool supports multiple sequencing platforms and provides detailed hit statistics
 * for each target gene, enabling rapid serotype identification for epidemiological investigations
 * and outbreak response.
 *
 * @status stable
 * @keywords Shigella, serotyping, k-mer, Illumina, Nanopore, epidemiology, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, shigatyper
 *
 * @subworkflows utils_bactopia-tools, shigatyper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.tsv                              Tab-delimited file containing predicted Shigella serotype for each sample
 * @publish *-hits.tsv                         Detailed statistics about each individual gene hit used for serotype prediction
 *
 * @section Merged Results
 * @publish shigatyper.tsv                     Merged TSV file containing ShigaTyper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                            Tool execution logs including ShigaTyper output
 * @publish logs/nf-*                          Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml                       Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { SHIGATYPER          } from '../../../subworkflows/shigatyper/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_shigatyper = SHIGATYPER(ch_bactopiatool.reads)

    publish:
    // Per-sample
    sample_outputs = ch_shigatyper.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_shigatyper.sample_outputs)
    // Run-level
    run_outputs = ch_shigatyper.run_outputs
    run_nf_logs = collectNextflowLogs(ch_shigatyper.run_outputs)
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
