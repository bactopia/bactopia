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
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    KLEBORATE(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(KLEBORATE.out.results)
    ch_logs = ch_logs.mix(KLEBORATE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(KLEBORATE.out.nf_logs)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions)

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
