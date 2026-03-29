#!/usr/bin/env nextflow
/**
 * Identify insertion sequence positions in bacterial genomes.
 *
 * This Bactopia Tool uses [ISMapper](https://github.com/jhawkey/IS_mapper) to identify
 * transposase insertion sites in bacterial genomes from short read sequence data.
 *
 * @status stable
 * @keywords insertion sequences, transposons, comparative genomics, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,comparative
 * @citation ismapper
 *
 * @subworkflows bactopiatool_init, ismapper
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input reference
 * Reference genome in FASTA format
 *
 * @input insertions
 * File containing insertion sequences to search for
 *
 * @section Per-Sample Results
 * @publish *.vcf            VCF file with insertion site calls
 * @publish *.txt            Summary of insertion site findings
 *
 * @section Merged Results
 * @publish ismapper.tsv      Merged TSV file containing ISMapper results from all samples
 *
 * @section Execution Logs
 * @publish logs/**           Tool execution logs
 * @publish logs/nf-*         Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml      Software version information
 */
nextflow.preview.types = true

params {
    rundir : String

    // Tool-specific parameters
    reference  : Path
    insertions : Path
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { ISMAPPER          } from '../../../subworkflows/ismapper/main'

include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    BACTOPIATOOL_INIT()
    ISMAPPER(
        BACTOPIATOOL_INIT.out.reads,
        params.reference,
        params.insertions
    )

    ch_sample_nf_logs = collectNextflowLogs(ISMAPPER.out.sample_outputs)
    ch_run_nf_logs = collectNextflowLogs(ISMAPPER.out.run_outputs)

    publish:
    // Per-sample records (scope: sample)
    sample_outputs = ISMAPPER.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    // Run-level records (scope: run)
    run_outputs = ISMAPPER.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_outputs {
        path { r ->
            r.supplemental >> "${r.meta.output_dir}/"
            r.logs         >> "${r.meta.logs_dir}/"
            r.versions     >> "${r.meta.logs_dir}/"
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
