//
// pangenome - Pangenome analysis with optional core-genome phylogeny
//

include { NCBIGENOMEDOWNLOAD } from '../ncbigenomedownload/main' addParams( options: [] )
//include { PROKKA } from '../../../modules/nf-core/modules/prokka/main' addParams( options: [] )

if (params.use_roary) {
    include { ROARY as PG_TOOL } from '../../../modules/nf-core/modules/roary/main' addParams( options: [] )
} else {
    include { PIRATE as PG_TOOL } from '../../../modules/nf-core/modules/pirate/main' addParams( options: [ args: "", is_module: false] )
}

include { IQTREE as START_TREE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: [] )
include { CLONALFRAMEML } from '../../../modules/nf-core/modules/clonalframeml/main' addParams( options: [] )
include { IQTREE as FINAL_TREE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: [] )
include { SNPDISTS } from '../../../modules/nf-core/modules/snpdists/main' addParams( options: [] )
include { SCOARY } from '../../../modules/nf-core/modules/scoary/main' addParams( options: [] )
                                                                       
workflow PANGENOME {
    take:
    gff // channel: [ val(meta), [ gff ] ]

    main:
    ch_versions = Channel.empty()
    ch_needs_prokka = Channel.empty()

    // Include public genomes (optional)
    /*
    if (params.accession || params.accessions || params.species) {
        NCBIGENOMEDOWNLOAD()
        NCBIGENOMEDOWNLOAD.out.fna.collect{fna -> fna}.map{ fna -> [[id:fna.getSimpleName()], fna]}.set{ ch_to_prokka }
        ch_needs_prokka.mix(ch_to_prokka)
        ch_versions.mix(NCBIGENOMEDOWNLOAD.out.versions.first())
    }
    */

    // Collect local assemblies
    /*
    if (params.assembly) {
        assemblies = []
        if (file(params.assembly).exists()) {
            if (file(params.assembly).isDirectory()) {
                assemblies_found = file("${params.assembly}/${params.assembly_pattern}")
                if (assemblies_found.size() == 0) {
                    log.error("0 assemblies were found in ${params.assembly} using the pattern ${params.assembly_pattern}, please check. Unable to continue.")
                    exit 1
                } else {
                    assemblies_found.each { assembly ->
                        assemblies << [[id: file(assembly).getSimpleName()], file(assembly).getSimpleName()]
                    }
                }
            } else {
                assemblies << [[id: file(params.assembly).getSimpleName()], file(params.assembly).getSimpleName()]
                is_compressed = params.assembly.endsWith(".gz") ? true : false
                has_assembly = true
            }
        } else {
            log.error("Could not open ${params.assembly}, please verify existence. Unable to continue.")
            exit 1
        }
        log.info("Found ${assemblies.size()} local assemblies.")
    }
    */

    // Annotate non-Bactopia genomes
    //if (ch_needs_prokka.count()) {
    //    PROKKA(ch_needs_prokka)
    ///    gff.mix(PROKKA.out.gff)
    //    ch_versions.mix(PROKKA.out.versions.first())
    //}

    // Create Pangenome
    PG_TOOL(gff)
    ch_versions.mix(PG_TOOL.out.versions.first())

    // Per-sample SNP distances
    SNPDISTS(PG_TOOL.out.aln)
    ch_versions.mix(PG_TOOL.out.versions.first())
    
    // Identify Recombination
    if (!params.skip_recombination) {
        START_TREE(PG_TOOL.out.aln)
        CLONALFRAMEML(START_TREE.out.tree, PG_TOOL.out.aln)
        ch_versions.mix(CLONALFRAMEML.out.versions.first())
    }

    // Create core-genome phylogeny
    if (!params.skip_phylogeny) {
        if (params.skip_recombination) {
            FINAL_TREE(PG_TOOL.out.aln)
        } else {
            FINAL_TREE(CLONALFRAMEML.out.masked_aln)
        }
        ch_versions.mix(FINAL_TREE.out.versions.first())
    }

    // Pan-genome GWAS
    if (params.traits) {
        SCOARY([id:'scoary'], PG_TOOL.out.csv, file(params.traits))
        ch_versions.mix(SCOARY.out.versions.first())
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
