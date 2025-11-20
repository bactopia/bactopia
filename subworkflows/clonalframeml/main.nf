//
// clonalframeml - Predict recombination events in bacterial genomes
//
include { CLONALFRAMEML as CLONALFRAMEML_MODULE } from '../../modules/clonalframeml/main'
include { IQTREE                                } from '../iqtree/main'
include { SNPDISTS                              } from '../snpdists/main'

workflow CLONALFRAMEML {
    take:
    alignment // channel: [ val(meta), [ aln ] ]

    main:
    ch_results = channel.empty()
    ch_logs = channel.empty()
    ch_nf_logs = channel.empty()
    ch_versions = channel.empty()
    
    // Create a quick start tree
    alignment.collect{_meta, aln -> aln}.map{ aln -> [[name: "iqtree-fast", process_name: "iqtree-fast"], aln]}.set{ ch_aln }
    IQTREE(ch_aln)
    ch_results = ch_results.mix(IQTREE.out.results)
    ch_logs = ch_logs.mix(IQTREE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(IQTREE.out.nf_logs)
    ch_versions = ch_versions.mix(IQTREE.out.versions)

    // Run ClonalFrameML
    IQTREE.out.aln_tree.collect{_meta, aln, treefile -> [aln, treefile]}.map{ aln, treefile -> [[name: "core-genome", process_name: "clonalframeml"], aln, treefile]}.set{ ch_aln_tree }
    CLONALFRAMEML_MODULE(ch_aln_tree)
    ch_results = ch_results.mix(
        CLONALFRAMEML_MODULE.out.emsim,
        CLONALFRAMEML_MODULE.out.em,
        CLONALFRAMEML_MODULE.out.status,
        CLONALFRAMEML_MODULE.out.newick,
        CLONALFRAMEML_MODULE.out.fasta,
        CLONALFRAMEML_MODULE.out.pos_ref,
        CLONALFRAMEML_MODULE.out.masked_aln
    )
    ch_versions = ch_versions.mix(CLONALFRAMEML_MODULE.out.versions)
    ch_logs = ch_logs.mix(CLONALFRAMEML_MODULE.out.logs)

    // Per-sample SNP distances
    CLONALFRAMEML_MODULE.out.masked_aln.collect{_meta, aln -> aln}.map{ aln -> [[name: "core-genome.masked.distance", process_name: "snpdists-masked"], aln]}.set{ ch_masked_aln }
    SNPDISTS(ch_masked_aln)
    ch_results = ch_results.mix(SNPDISTS.out.tsv)
    ch_logs = ch_logs.mix(SNPDISTS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNPDISTS.out.nf_logs)
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)

    emit:
    // Individual outputs
    masked_aln = CLONALFRAMEML_MODULE.out.masked_aln

    // Generic aggregate outputs
    results = ch_results
    logs = ch_logs
    nf_logs = ch_nf_logs.mix(
        CLONALFRAMEML_MODULE.out.nf_begin,
        CLONALFRAMEML_MODULE.out.nf_err,
        CLONALFRAMEML_MODULE.out.nf_log,
        CLONALFRAMEML_MODULE.out.nf_out,
        CLONALFRAMEML_MODULE.out.nf_run,
        CLONALFRAMEML_MODULE.out.nf_sh,
        CLONALFRAMEML_MODULE.out.nf_trace
    )
    versions = ch_versions
}
