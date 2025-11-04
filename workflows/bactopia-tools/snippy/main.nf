#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { SNIPPY             } from '../../../subworkflows/snippy/run/main'
include { SNIPPY_CORE        } from '../../../subworkflows/snippy/core/main'
include { GUBBINS            } from '../../../subworkflows/gubbins/main'
include { IQTREE             } from '../../../subworkflows/iqtree/main'
include { paramsHelp         } from 'plugin/nf-bactopia'
include { workflowSummary    } from 'plugin/nf-bactopia'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow {

    main:
    // Check if help is requested
    if (params.help || params.help_all) {
        log.info paramsHelp()
        exit 0
    }

    // Initialize and execute the workflow
    ch_results = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()
    ch_versions = Channel.empty()
    BACTOPIATOOL_INIT(params.bactopia, params.workflow.ext, params.include, params.exclude)

    // Download if applicable
    ch_reference = Channel.empty()
    if (params.reference) {
        ch_reference = [[id: 'snippy'], file(params.reference)]
    } else if (params.accession) {
        NCBIGENOMEDOWNLOAD([])
        ch_reference = NCBIGENOMEDOWNLOAD.out.bactopia_tools.first()
    }

    // Run Snippy per-sample
    SNIPPY(BACTOPIATOOL_INIT.out.samples, ch_reference)
    ch_results = ch_results.mix(SNIPPY.out.results)
    ch_logs = ch_logs.mix(SNIPPY.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNIPPY.out.nf_logs)
    ch_versions = ch_versions.mix(SNIPPY.out.versions)

    // Identify core SNPs
    SNIPPY.out.vcf.collect{_meta, vcf -> vcf}.map{ vcf -> [[id:'core-snp'], vcf]}.set{ ch_merge_vcf }
    SNIPPY.out.aligned_fa.collect{_meta, aligned_fa -> aligned_fa}.map{ aligned_fa -> [[id:'core-snp'], aligned_fa]}.set{ ch_merge_aligned_fa }
    ch_merge_vcf.join( ch_merge_aligned_fa ).set{ ch_snippy_core }

    // Identify core SNPs
    SNIPPY_CORE(ch_snippy_core, ch_reference, params.snippy_core_mask ? file(params.snippy_core_mask, checkIfExists: true) : [])
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

    workflow.onComplete {
        log.info workflowSummary()
    }

    publish:
    results = ch_results
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions
}

output {
    results {
        path { meta, _file -> "${meta.output_dir}/" }
    }
    logs {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
    nf_logs {
        path { meta, file -> {
            file >> "${meta.logs_dir}/nf${file.name}"
        } }
    }
    versions {
        path { meta, _file -> "${meta.logs_dir}/" }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
