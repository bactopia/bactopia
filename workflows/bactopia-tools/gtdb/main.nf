#!/usr/bin/env nextflow
/**
 * Identify marker genes and assign taxonomic classifications using GTDB.
 *
 * This Bactopia Tool uses [GTDB-Tk's](https://github.com/Ecogenomics/GTDBTk) classify
 * workflow to assign taxonomic classifications to samples using the
 * [Genome Taxonomy Database](https://gtdb.ecogenomic.org/).
 *
 * @status stable
 * @keywords taxonomy, classification, marker genes, phylogeny, gtdb, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,taxonomy,database-dependent
 * @citation gtdb
 *
 * @subworkflows bactopiatool_init, gtdb
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input gtdb
 * Path to GTDB database directory
 *
 * @input download_gtdb
 * Download GTDB database if not found
 *
 * @input gtdb_save_as_tarball
 * Save GTDB database as tarball after download
 *
 * @section Taxonomic Classification
 * @publish *.summary.tsv        Taxonomic classification summary
 * @publish *.gtdbtk.tsv        Detailed GTDB-Tk results
 *
 * @section Merged Results
 * @publish gtdb.tsv            Merged TSV file containing GTDB results from all samples
 *
 * @section Execution Logs
 * @publish logs/**             Tool execution logs
 * @publish logs/nf-*           Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml        Software version information
 */
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    gtdb                 : Path
    download_gtdb        : Boolean
    gtdb_save_as_tarball : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { GTDB              } from '../../../subworkflows/gtdb/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Set<Path>>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Set<Path>>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    GTDB(
        BACTOPIATOOL_INIT.out.assembly,
        params.gtdb,
        params.download_gtdb,
        params.gtdb_save_as_tarball
    )

    // Collect outputs
    ch_results = ch_results.mix(GTDB.out.results)
    ch_logs = ch_logs.mix(GTDB.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GTDB.out.nf_logs)
    ch_versions = ch_versions.mix(GTDB.out.versions)

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
