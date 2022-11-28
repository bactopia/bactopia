//
// gubbins - Rapid phylogenetic analysis of recombinant bacterial sequences
//
include { initOptions } from '../../../lib/nf/functions'
options = initOptions(params.containsKey("options") ? params.options : [:], 'gubbins')
options.is_module = params.wf == 'gubbins' ? true : false
options.suffix = options.suffix ? options.suffix : 'gubbins'
options.args = [
    params.remove_identical_sequences ? "--remove-identical-sequences" : "",
    "--iterations ${params.iterations}",
    "--min-snps ${params.min_snps}",
    "--min-window-size ${params.min_window_size}",
    "--max-window-size ${params.max_window_size}",
    "--filter-percentage ${params.filter_percentage }",
    "${params.gubbin_opts}"
].join(' ').replaceAll("\\s{2,}", " ").trim()

include { GUBBINS as GUBBINS_MODULE } from '../../../modules/nf-core/gubbins/main' addParams( options: options )

workflow GUBBINS {
    take:
    alignment // channel: [ val(meta), [ aln ] ]

    main:
    ch_versions = Channel.empty()

    GUBBINS_MODULE(alignment)
    ch_versions = ch_versions.mix(GUBBINS_MODULE.out.versions)

    emit:
    masked_aln = GUBBINS_MODULE.out.masked_aln
    fasta = GUBBINS_MODULE.out.fasta
    gff = GUBBINS_MODULE.out.gff
    vcf = GUBBINS_MODULE.out.vcf
    phylip = GUBBINS_MODULE.out.phylip
    embl_predicted = GUBBINS_MODULE.out.embl_predicted
    embl_branch = GUBBINS_MODULE.out.embl_branch
    tree = GUBBINS_MODULE.out.tree
    tree_labelled = GUBBINS_MODULE.out.tree_labelled
    versions = ch_versions // channel: [ versions.yml ]
}
