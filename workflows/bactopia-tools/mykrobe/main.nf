#!/usr/bin/env nextflow
/**
 * Antimicrobial resistance detection for specific bacterial species.
 *
 * This Bactopia Tool uses [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) to predict
 * antimicrobial resistance for _Mycobacterium tuberculosis_, _Staphylococcus aureus_,
 * _Shigella sonnei_, and _Salmonella typhi_ from sequencing data.
 *
 * @status stable
 * @keywords fastq, antimicrobial resistance, species-specific, mykrobe, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,species-specific
 * @citation mykrobe
 *
 * @subworkflows bactopiatool_init, mykrobe
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input mykrobe_species
 * Target species for AMR prediction (e.g., tb, saureus, sonnei, typhi)
 *
 * @section Per-Sample Results
 * @publish *.json            Mykrobe analysis results in JSON format
 * @publish *.txt             Tab-delimited report of resistance predictions
 *
 * @section Merged Results
 * @publish mykrobe.tsv       Merged TSV file containing results from all samples
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
    rundir   : String

    // Tool-specific parameters
    mykrobe_species : String
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { MYKROBE           } from '../../../subworkflows/mykrobe/main'

workflow {
    main:
    BACTOPIATOOL_INIT()
    MYKROBE(
        BACTOPIATOOL_INIT.out.reads,
        params.mykrobe_species
    )
    ch_sample_nf_logs = MYKROBE.out.sample_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }
    ch_run_nf_logs = MYKROBE.out.run_outputs.flatMap { r ->
        r.nf_logs.collect { f -> tuple(r.meta, f) }
    }

    publish:
    sample_outputs = MYKROBE.out.sample_outputs
    sample_nf_logs = ch_sample_nf_logs
    run_outputs = MYKROBE.out.run_outputs
    run_nf_logs = ch_run_nf_logs
}

output {
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
