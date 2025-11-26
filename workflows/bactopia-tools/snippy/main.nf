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
    reference          : Path
    accession          : String
    snippy_core_mask   : Path
    skip_recombination : Boolean
    skip_phylogeny     : Boolean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools/main'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { SNIPPY             } from '../../../subworkflows/snippy/run/main'
include { SNIPPY_CORE        } from '../../../subworkflows/snippy/core/main'
include { GUBBINS            } from '../../../subworkflows/gubbins/main'
include { IQTREE             } from '../../../subworkflows/iqtree/main'
include { formatSamples      } from 'plugin/nf-bactopia'

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

    // Download if applicable
    ch_reference = channel.empty()
    if (params.reference) {
        ch_reference = [[id: 'snippy'], params.reference]
    } else if (params.accession) {
        NCBIGENOMEDOWNLOAD([])
        ch_reference = NCBIGENOMEDOWNLOAD.out.bactopia_tools.first()
    }

    // Run Snippy per-sample
    SNIPPY(
        formatSamples(BACTOPIATOOL_INIT.out.samples, BACTOPIATOOL_INIT.out.data_types),
        ch_reference
    )
    ch_results = ch_results.mix(SNIPPY.out.results)
    ch_logs = ch_logs.mix(SNIPPY.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNIPPY.out.nf_logs)
    ch_versions = ch_versions.mix(SNIPPY.out.versions)

    // Identify core SNPs
    ch_merge_vcf = SNIPPY.out.vcf.collect{_meta, vcf -> vcf}.map{ vcf -> [[id:'core-snp'], vcf]}
    ch_merge_aligned_fa = SNIPPY.out.aligned_fa.collect{_meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'core-snp'], aligned_fa]}
    ch_snippy_core = ch_merge_vcf.join( ch_merge_aligned_fa )

    // Identify core SNPs
    SNIPPY_CORE(ch_snippy_core, ch_reference, params.snippy_core_mask ? [params.snippy_core_mask] : [])
    ch_results = ch_results.mix(SNIPPY_CORE.out.results)
    ch_logs = ch_logs.mix(SNIPPY_CORE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNIPPY_CORE.out.nf_logs)
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)

    // (optional) Identify Recombination
    if (!params.skip_recombination) {
        // Run Gubbins
        GUBBINS(SNIPPY_CORE.out.clean_full_aln)
        ch_results = ch_results.mix(GUBBINS.out.results)
        ch_logs = ch_logs.mix(GUBBINS.out.logs)
        ch_nf_logs = ch_nf_logs.mix(GUBBINS.out.nf_logs)
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
    }

    // Create core-snp phylogeny
    if (!params.skip_phylogeny) {
        if (!params.skip_recombination) {
            IQTREE(GUBBINS.out.masked_aln)
        } else {
            IQTREE(SNIPPY_CORE.out.clean_full_aln)
        }
        ch_results = ch_results.mix(IQTREE.out.results)
        ch_logs = ch_logs.mix(IQTREE.out.logs)
        ch_nf_logs = ch_nf_logs.mix(IQTREE.out.nf_logs)
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

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
