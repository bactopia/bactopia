/**
 * Detect and mask recombination events in bacterial phylogenies.
 *
 * This subworkflow uses [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML) to
 * detect and mask recombination events in bacterial phylogenies. It first builds a quick
 * phylogenetic tree using [IQ-TREE](https://github.com/iqtree/iqtree2), then identifies
 * recombination regions and creates a recombination-masked alignment. Finally, it
 * calculates SNP distances from the masked alignment using [snp-dists](https://github.com/tseemann/snp-dists).
 *
 * @status stable
 * @keywords recombination, phylogeny, masking, clonalframe, bacterial evolution
 * @tags complexity:moderate input-type:single output-type:multiple features:aggregation
 * @citation clonalframeml, iqtree, snpdists
 *
 * @subworkflows iqtree, snpdists
 * @modules clonalframeml
 *
 * @input tuple(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Core-genome alignment in FASTA format
 *
 * @output masked_aln     Recombination-masked alignment file with detected recombination regions removed
 * @output results        Aggregated results channel containing all output files from ClonalFrameML, IQ-TREE, and snp-dists
 * @output logs           Aggregated logs channel containing all execution logs
 * @output nf_logs        Aggregated Nextflow execution scripts and logs for debugging from all processes
 * @output versions       Aggregated version information from all executed tools
 */
nextflow.preview.types = true

include { CLONALFRAMEML as CLONALFRAMEML_MODULE } from '../../modules/clonalframeml/main'
include { IQTREE                                } from '../iqtree/main'
include { SNPDISTS                              } from '../snpdists/main'
include { flattenPaths                          } from 'plugin/nf-bactopia'
include { gather                                } from 'plugin/nf-bactopia'

workflow CLONALFRAMEML {
    take:
    alignment: Channel<Tuple<Map, Set<Path>>>

    main:
    ch_results = channel.empty() as Channel<Tuple<Map, Path>>
    ch_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_nf_logs = channel.empty() as Channel<Tuple<Map, Path>>
    ch_versions = channel.empty() as Channel<Tuple<Map, Path>>

    // Execute subworkflows
    // Create a quick start tree
    IQTREE(gather(alignment, 'iqtree-fast', 'name: "iqtree-fast", process_name: "iqtree-fast"'))
    ch_results = ch_results.mix(IQTREE.out.results)
    ch_logs = ch_logs.mix(IQTREE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(IQTREE.out.nf_logs)
    ch_versions = ch_versions.mix(IQTREE.out.versions)

    // Run ClonalFrameML - gather alignment and tree together
    ch_aln_tree = IQTREE.out.aln_tree.collect{_meta, aln, treefile -> [aln, treefile]}.map{ aln, treefile -> [[name: "core-genome", process_name: "clonalframeml"], aln, treefile]}
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
    SNPDISTS(gather(CLONALFRAMEML_MODULE.out.masked_aln, 'core-genome.masked.distance', 'name: "core-genome.masked.distance", process_name: "snpdists-masked"'))
    ch_results = ch_results.mix(SNPDISTS.out.tsv)
    ch_logs = ch_logs.mix(SNPDISTS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SNPDISTS.out.nf_logs)
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)

    emit:
    // Individual outputs
    masked_aln: Channel<Tuple<Map, Path>> = CLONALFRAMEML_MODULE.out.masked_aln

    // Generic aggregate outputs
    results: Channel<Tuple<Map, Path>> = flattenPaths([ch_results])
    logs: Channel<Tuple<Map, Path>> = flattenPaths([ch_logs])
    nf_logs: Channel<Tuple<Map, Path>> = flattenPaths([
        ch_nf_logs,
        CLONALFRAMEML_MODULE.out.nf_logs
    ])
    versions: Channel<Tuple<Map, Path>> = flattenPaths([ch_versions])
}
