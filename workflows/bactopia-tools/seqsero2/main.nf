#!/usr/bin/env nextflow
/**
 * Salmonella serotype prediction from sequencing reads or assemblies.
 *
 * This Bactopia Tool uses [SeqSero2](https://github.com/denglab/SeqSero2) to predict Salmonella
 * serotypes from both raw sequencing reads and assembled genomes. SeqSero2 is a novel pipeline
 * for determining Salmonella serotypes using raw sequencing reads or assemblies through
 * k-mer analysis and targeted identification of O and H antigen genes. The tool provides
 * accurate serotype predictions following the Kaufmann-White scheme, supporting traditional
 * and molecular serotyping methods for epidemiological surveillance and outbreak investigation.
 *
 * @status stable
 * @keywords Salmonella, serotyping, epidemiology, O antigen, H antigen, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,aggregation
 * @citation csvtk, seqsero2
 *
 * @subworkflows bactopiatool_init, seqsero2
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *_result.tsv             Tab-delimited file with detailed SeqSero2 results for each sample
 * @publish *_result.txt             Text file with key-value pairs of SeqSero2 prediction results
 * @publish *_log.txt                Detailed log file from SeqSero2 analysis
 *
 * @section Merged Results
 * @publish seqsero2.tsv             Merged TSV file containing SeqSero2 results from all samples
 *
 * @section Execution Logs
 * @publish logs/**                  Tool execution logs including SeqSero2 output
 * @publish logs/nf-*                Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml             Software version information
 */
nextflow.preview.types = true

params {
    rundir : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { SEQSERO2          } from '../../../subworkflows/seqsero2/main'

include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    SEQSERO2(BACTOPIATOOL_INIT.out.assembly)

    ch_sample_nf_logs = collectNextflowLogs(SEQSERO2.out.sample_outputs)
    ch_run_nf_logs = collectNextflowLogs(SEQSERO2.out.run_outputs)

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = SEQSERO2.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = SEQSERO2.out.run_outputs
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
