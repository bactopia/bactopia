#!/usr/bin/env nextflow
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BACTOPIATOOL_INIT  } from '../../../subworkflows/utils/bactopia-tools'
include { NCBIGENOMEDOWNLOAD } from '../../../subworkflows/ncbigenomedownload/main'
include { PROKKA             } from '../../../subworkflows/prokka/main'
include { PANGENOME          } from '../../../subworkflows/pangenome/main'
include { CLONALFRAMEML      } from '../../../subworkflows/clonalframeml/main'
include { IQTREE             } from '../../../subworkflows/iqtree/main'
include { SCOARY             } from '../../../subworkflows/scoary/main'
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
    ch_samples = BACTOPIATOOL_INIT.out.samples

    // Download if applicable
    if (params.species || params.accession || params.accessions) {
        NCBIGENOMEDOWNLOAD(params.accessions ? file(params.accessions) : [])
        PROKKA(
            NCBIGENOMEDOWNLOAD.out.bactopia_tools,
            params.proteins ? file(params.proteins) : [],
            params.prodigal_tf ? file(params.prodigal_tf) : []
        )
        ch_samples = ch_samples.mix(PROKKA.out.gff)
    }

    // Create the pangenome
    ch_samples.collect{_meta, gff -> gff}.map{ gff -> [[id: params.use_pirate? 'pirate' : (params.use_roary ? 'roary' : 'panaroo')], gff]}.set{ ch_merge_gff }
    PANGENOME(ch_merge_gff, params.use_pirate, params.use_roary)
    ch_results = ch_results.mix(PANGENOME.out.results)
    ch_logs = ch_logs.mix(PANGENOME.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PANGENOME.out.nf_logs)
    ch_versions = ch_versions.mix(PANGENOME.out.versions)
    
    // (optional) Identify Recombination
    if (!params.skip_recombination) {
        // Run ClonalFrameML
        CLONALFRAMEML(PANGENOME.out.aln)
        ch_results = ch_results.mix(CLONALFRAMEML.out.results)
        ch_logs = ch_logs.mix(CLONALFRAMEML.out.logs)
        ch_nf_logs = ch_nf_logs.mix(CLONALFRAMEML.out.nf_logs)
        ch_versions = ch_versions.mix(CLONALFRAMEML.out.versions)
    }

    // (optional) Create core-genome phylogeny
    if (!params.skip_phylogeny) {
        ch_final_aln = Channel.empty()
        if (params.skip_recombination) {
            PANGENOME.out.aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-genome", process_name: "iqtree"], aln]}.set{ ch_final_aln }
        } else {
            CLONALFRAMEML.out.masked_aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-genome", process_name: "iqtree"], aln]}.set{ ch_final_aln }
        }
        IQTREE(ch_final_aln)
        ch_results = ch_results.mix(IQTREE.out.results)
        ch_logs = ch_logs.mix(IQTREE.out.logs)
        ch_nf_logs = ch_nf_logs.mix(IQTREE.out.nf_logs)
        ch_versions = ch_versions.mix(IQTREE.out.versions)
    }

    // Pan-genome GWAS
    if (params.traits) {
        SCOARY(PANGENOME.out.csv, file(params.traits))
        ch_results = ch_results.mix(SCOARY.out.csv)
        ch_logs = ch_logs.mix(SCOARY.out.logs)
        ch_nf_logs = ch_nf_logs.mix(SCOARY.out.nf_logs)
        ch_versions = ch_versions.mix(SCOARY.out.versions)
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
