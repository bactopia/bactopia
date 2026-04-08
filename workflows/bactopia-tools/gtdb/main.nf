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
    rundir : String

    // Tool-specific parameters
    gtdb                 : Value<Path>
    download_gtdb        : Boolean
    gtdb_save_as_tarball : Boolean
}

include { BACTOPIATOOL_INIT   } from '../../../subworkflows/utils/bactopia-tools/main'
include { GTDB                } from '../../../subworkflows/gtdb/main'
include { collectNextflowLogs } from 'plugin/nf-bactopia'

workflow {
    main:
    ch_bactopiatool = BACTOPIATOOL_INIT()
    ch_gtdb = GTDB(
        ch_bactopiatool.assembly,
        params.gtdb,
        params.download_gtdb,
        params.gtdb_save_as_tarball
    )

    publish:
    // Per-sample
    sample_outputs = ch_gtdb.sample_outputs
    sample_nf_logs = collectNextflowLogs(ch_gtdb.sample_outputs)
    // Run-level
    run_outputs = ch_gtdb.run_outputs
    run_nf_logs = collectNextflowLogs(ch_gtdb.run_outputs)
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
