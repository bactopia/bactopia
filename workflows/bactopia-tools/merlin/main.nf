#!/usr/bin/env nextflow
nextflow.preview.types = true
/**
 * Bactopia Tool: Merlin.
 *
 * MinmER assisted species-specific bactopia tool seLectIoN
 * _MinmER assisted species-specific bactopia tool seLectIoN_, or Merlin, uses distances based
 * on the RefSeq sketch downloaded by `bactopia datasets` to automatically run species-specific tools.
 * Currently Merlin knows 16 spells for which cover the following:
 * | Genus/Species | Tools |
 * |---------------|-------|
 * | Escherichia / Shigella   | [ECTyper](../bactopia-tools/ectyper.md), [ShigaTyper](../bactopia-tools/shigatyper.md), [ShigEiFinder](../bactopia-tools/shigeifinder.md)  |
 * | Haemophilus   | [hicap](../bactopia-tools/hicap.md), [HpsuisSero](../bactopia-tools/ssuissero.md) |
 * | Klebsiella | [Kleborate](../bactopia-tools/kleborate.md) |
 * | Legionella | [legsta](../bactopia-tools/legsta.md) |
 * | Listeria | [LisSero](../bactopia-tools/lissero.md) |
 * | Mycobacterium | [TBProfiler](../bactopia-tools/tbprofiler.md) |
 * | Neisseria | [meningotype](../bactopia-tools/meningotype.md), [ngmaster](../bactopia-tools/ngmaster.md) |
 * | Pseudomonas | [pasty](../bactopia-tools/pasty.md) |
 * | Salmonella | [SeqSero2](../bactopia-tools/seqsero2.md), [SISTR](../bactopia-tools/sistr.md) |
 * | Staphylococcus | [AgrVATE](../bactopia-tools/agrvate.md), [spaTyper](../bactopia-tools/spatyper.md), [staphopia-sccmec](../bactopia-tools/staphopiasccmec.md) |
 * | Streptococcus | [emmtyper](../bactopia-tools/emmtyper.md), [pbptyper](../bactopia-tools/pbptyper.md), [SsuisSero](../bactopia-tools/ssuissero.md) |
 * Merlin is avialable as an independent Bactopia Tool, or in the Bactopia with the `--ask_merlin` parameter. Even better,
 * if you want to force Merlin to execute all species-specific tools (no matter the distance), you can use `--full_merlin`.
 * Then all the spells will be unleashed!
 *
 * @status stable
 * @keywords serotype, species-specific
 *
 * @subworkflows bactopiatool_init, merlin
 *
 * @input rundir
 * Run directory containing Bactopia results
 *
 * @section Per-Sample Results
 * @publish *    Analysis results
 *
 * @section Merged Results
 * @publish merged-*    Aggregated results from all samples
 *
 * @section Execution Logs
 * @publish logs/**   Tool execution logs
 * @publish logs/nf-* Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml Software version information
   */

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    emmtyper_blastdb      : Path?
    hicap_database_dir    : Path?
    hicap_model_fp        : Path?
    spatyper_repeats      : Path?
    spatyper_repeat_order : Path?
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { DATASETS          } from '../../../modules/bactopia/datasets/main'
include { MERLIN            } from '../../../subworkflows/merlin/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    DATASETS()
    MERLIN(
        BACTOPIATOOL_INIT.out.samples_2,
        DATASETS.out.mash_db,
        // emmtyper
        params.emmtyper_blastdb,
        // hicap
        params.hicap_database_dir,
        params.hicap_model_fp,
        // staphtyper
        params.spatyper_repeats,
        params.spatyper_repeat_order
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
    run_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.output_dir}" }
    }
    run_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }
    run_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${params.rundir}/${meta.logs_dir}/nf${file.name}"
        }
    }
    run_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${params.rundir}/${meta.logs_dir}/" }
    }

    // Sample-level outputs (stored in ${params.outdir}/<SAMPLE_NAME>/)
    sample_results: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    sample_logs: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    sample_nf_logs: Channel<Tuple<Map, Path>> {
        path { meta, file ->
            file >> "${meta.logs_dir}/nf${file.name}"
        }
    }
    sample_versions: Channel<Tuple<Map, Path>> {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}
