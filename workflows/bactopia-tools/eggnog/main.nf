#!/usr/bin/env nextflow
/**
 * Functional annotation of proteins using orthologous groups and phylogenies.
 *
 * This Bactopia Tool uses [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) to assign
 * functional annotation to protein sequences. eggNOG-mapper uses orthologous groups and phylogenies
 * from the eggNOG database to more precisely functionally annotate than traditional homology methods.
 *
 * @status stable
 * @keywords functional annotation, orthology, proteins, eggnog, bactopia-tool
 * @tags complexity:complex input-type:parameter output-type:multiple features:bactopia-tool,functional-annotation,orthology,database-dependent
 * @citation eggnog
 *
 * @subworkflows bactopiatool_init, eggnog
 *
 * @input rundir
 * Directory containing results from a completed Bactopia analysis run
 *
 * @input eggnog_db
 * Path to eggNOG database for functional annotation
 *
 * @input download_eggnog
 * Download eggNOG database if not found locally
 *
 * @section Annotation
 * @publish *.emapper.annotations      Results from the annotation phase
 * @publish *.emapper.hits             Results from the search phase (HMMER, Diamond or MMseqs2)
 * @publish *.emapper.seed_orthologs   Results from parsing the hits
 * @publish *.emapper.annotations.xlsx Annotations in Excel format
 * @publish *.emapper.orthologs        List of orthologs found for each query
 * @publish *.emapper.genepred.fasta   Sequences of predicted CDS
 * @publish *.emapper.gff              GFF of predicted CDS
 * @publish *.emapper.no_annotations.fasta Sequences without annotation
 * @publish *.emapper.pfam             Positions of PFAM domains identified
 *
 * @section Merged Results
 * @publish eggnog.tsv                 Merged annotations from all samples
 *
 * @section Execution Logs
 * @publish logs/**                    Tool execution logs
 * @publish logs/nf-*                  Nextflow execution scripts and logs for debugging
 *
 * @section Versions
 * @publish versions.yml               Software version information
 */
nextflow.preview.types = true

params {
    bactopia : String
    includes : String
    excludes : String
    workflow : Map
    rundir   : String

    // Tool-specific parameters
    eggnog_db       : Path
    download_eggnog : Boolean
}

include { BACTOPIATOOL_INIT } from '../../../subworkflows/utils/bactopia-tools/main'
include { EGGNOG            } from '../../../subworkflows/eggnog/main'

workflow {
    main:
    // Initialize output channels
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    BACTOPIATOOL_INIT()
    EGGNOG(
        BACTOPIATOOL_INIT.out.samples,
        params.eggnog_db,
        params.download_eggnog
    )

    // Collect outputs
    ch_results = ch_results.mix(EGGNOG.out.results)
    ch_logs = ch_logs.mix(EGGNOG.out.logs)
    ch_nf_logs = ch_nf_logs.mix(EGGNOG.out.nf_logs)
    ch_versions = ch_versions.mix(EGGNOG.out.versions)

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
