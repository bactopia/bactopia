#!/usr/bin/env nextflow
/**
 * Comprehensive typing of Neisseria meningitidis.
 *
 * This Bactopia Tool uses [meningotype](https://github.com/MDU-PHL/meningotype)
 * for _in silico_ typing of _Neisseria meningitidis_ genomes. It uses assembly contigs
 * to determine the serotype, MLST, finetyping (_porA_, _fetA_, _porB_), and
 * Bexsero antigen sequence typing (BAST) (_fHbp_, _NHBA_, _NadA_, _PorA_).
 *
 * @status stable
 * @keywords neisseria meningitidis, serotyping, mlst, finetyping, fasta, bactopia-tool
 * @tags complexity:moderate input-type:parameter output-type:multiple features:bactopia-tool,typing,mlst
 * @citation meningotype
 *
 * @subworkflows bactopiatool_init, meningotype
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @section Per-Sample Results
 * @publish *.txt                Comprehensive typing report
 * @publish *-allele.tsv         Allele typing results
 * @publish *-mlst.tsv           MLST typing results
 *
 * @section Merged Results
 * @publish meningotype.tsv       Merged TSV file containing meningotype results from all samples
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
include { MENINGOTYPE       } from '../../../subworkflows/meningotype/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    MENINGOTYPE(BACTOPIATOOL_INIT.out.samples)

    // Collect outputs
    ch_results = ch_results.mix(MENINGOTYPE.out.results)
    ch_logs = ch_logs.mix(MENINGOTYPE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MENINGOTYPE.out.nf_logs)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions)

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
