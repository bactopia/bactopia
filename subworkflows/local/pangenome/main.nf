//
// pangenome - Pangenome analysis with optional core-genome phylogeny
//
include { NCBIGENOMEDOWNLOAD } from '../ncbigenomedownload/main'
//include { PROKKA } from '../../../modules/nf-core/modules/prokka/main' addParams( options: [] )

if (params.use_roary) {
    include { ROARY as PG_TOOL } from '../roary/main' addParams( options: [] )
} else {
    include { PIRATE as PG_TOOL } from '../pirate/main' addParams( options: [publish_to_base: [".aln.gz"]] )
}

include { IQTREE as START_TREE } from '../../../modules/nf-core/modules/iqtree/main' addParams( options: [args: "-m ${params.m} -fast", suffix: 'start-tree', process_name: 'clonalframeml'])
include { CLONALFRAMEML } from '../clonalframeml/main' addParams( options: [publish_to_base: [".masked.aln.gz"]] )
include { IQTREE as FINAL_TREE } from '../iqtree/main' addParams( options: [suffix: 'core-genome', publish_to_base: [".iqtree"]] )
include { SNPDISTS } from '../../../modules/nf-core/modules/snpdists/main' addParams( options: [suffix: 'core-genome.distance', publish_to_base: true] )
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
    gff.collect{meta, gff -> gff}.map{ gff -> [[id: params.use_roary ? 'roary' : 'pirate'], gff]}.set{ ch_merge_gff }
    PG_TOOL(ch_merge_gff)
    ch_versions.mix(PG_TOOL.out.versions)

    // Per-sample SNP distances
    SNPDISTS(PG_TOOL.out.aln)
    ch_versions.mix(SNPDISTS.out.versions)
    
    // Identify Recombination
    if (!params.skip_recombination) {
        // Create quick tree for ClonalFrameML
        PG_TOOL.out.aln.collect{meta, aln -> aln}.map{ aln -> [[id: 'start-tree'], aln]}.set{ ch_start_tree }
        START_TREE(ch_start_tree)

        // Run ClonalFrameML
        PG_TOOL.out.aln.collect{meta, aln -> aln}.map{ aln -> [ aln ]}.set{ ch_aln }
        START_TREE.out.phylogeny.collect{meta, phylogeny -> phylogeny}.map{ phylogeny -> [ phylogeny ]}.set{ ch_phylogeny }
        CLONALFRAMEML(Channel.fromList([[id:'clonalframeml']]).combine(ch_phylogeny).combine(ch_aln))
        ch_versions.mix(CLONALFRAMEML.out.versions)
    }

    // Create core-genome phylogeny
    if (!params.skip_phylogeny) {
        if (params.skip_recombination) {
            FINAL_TREE(PG_TOOL.out.aln)
        } else {
            FINAL_TREE(CLONALFRAMEML.out.masked_aln)
        }
        ch_versions.mix(FINAL_TREE.out.versions)
    }

    // Pan-genome GWAS
    if (params.traits) {
        SCOARY([id:'scoary'], PG_TOOL.out.csv, file(params.traits))
        ch_versions.mix(SCOARY.out.versions)
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
