#!/usr/bin/env nextflow
nextflow.preview.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW PARAMETERS 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    emmtyper_blastdb      : Path
    hicap_database_dir    : Path
    hicap_model_fp        : Path
    spatyper_repeats      : Path
    spatyper_repeat_order : Path
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { formatSamples     } from '../../../subworkflows/utils/generic/main'
include { DATASETS          } from '../../../modules/bactopia/datasets/main'
include { MERLIN            } from '../../../subworkflows/merlin/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    // Initialize and execute the workflow
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    BACTOPIATOOL_INIT()
    DATASETS()
    MERLIN(
        formatSamples(BACTOPIATOOL_INIT.out.samples, BACTOPIATOOL_INIT.out.data_types),
        DATASETS.out.mash_db,
        // emmtyper
        params.emmtyper_blastdb ? [params.emmtyper_blastdb] : [],
        // hicap
        params.hicap_database_dir ? [params.hicap_database_dir] : [],
        params.hicap_model_fp ? [params.hicap_model_fp] : [],
        // staphtyper
        params.spatyper_repeats ? [params.spatyper_repeats] : [],
        params.spatyper_repeat_order ? [params.spatyper_repeat_order] : []
    )

    // Collect outputs
    ch_results = ch_results.mix(MERLIN.out.results)
    ch_logs = ch_logs.mix(MERLIN.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MERLIN.out.nf_logs)
    ch_versions = ch_versions.mix(MERLIN.out.versions)

    // Branch the based on scope (sample or run)
    ch_final_results = ch_results.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_logs = ch_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_nf_logs = ch_nf_logs.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    ch_final_versions = ch_versions.branch{ meta, _file ->
        run: meta.scope == 'run'
        sample: meta.scope == 'sample'
    }

    publish:
    run_results = ch_final_results.run
    run_logs = ch_final_logs.run
    run_nf_logs = ch_final_nf_logs.run
    run_versions = ch_final_versions.run
    sample_results = ch_final_results.sample
    sample_logs = ch_final_logs.sample
    sample_nf_logs = ch_final_nf_logs.sample
    sample_versions = ch_final_versions.sample
}

output {
    // Run-level outputs (stored in ${params.outdir}/bactopia-runs/<RUN_NAME>/)
    run_results {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs {
        path { meta, file -> {
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        } }
    }
    run_versions {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    sample_versions {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
