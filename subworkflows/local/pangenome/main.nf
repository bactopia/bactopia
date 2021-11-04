//
// pangenome - Pangenome analysis with optional core-genome phylogeny
//

params.OPTS = [:]

include { CLONALFRAMEML } from '../../../modules/nf-core/modules/pirate/main' addParams( options: [] )
include { IQTREE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: [] )
include { NCBIGENOMEDOWNLOAD } from '../../../modules/nf-core/modules/ncbigenomedownload/main' addParams( options: [] )
include { PIRATE } from '../../../modules/nf-core/modules/pirate/main' addParams( options: [] )
include { PROKKA } from '../../../modules/nf-core/modules/prokka/main' addParams( options: [] )
include { ROARY } from '../../../modules/nf-core/modules/roary/main' addParams( options: [] )
include { SNPDISTS } from '../../../modules/nf-core/modules/snpdists/main' addParams( options: [] )


workflow PANGENOME {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE DESCRIPTION
    //
    MODULE ( INPUTS )
    ch_versions = ch_versions.mix(MODULE.out.versions.first())

    // Include public genomes (optional)

        // Download genomes
        // Annotate genomes

    // If applicable add GFFs to full set


    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
