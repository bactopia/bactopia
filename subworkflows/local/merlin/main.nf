//
// merlin - MinmER assisted species-specific bactopia tool seLectIoN
//

include { MERLINDIST } from '../mashdist/main';
include { ECTYPER } from '../ectyper/main';
include { EMMTYPER } from '../emmtyper/main';
include { HICAP } from '../hicap/main';
include { HPSUISSERO } from '../hpsuissero/main';
include { GENOTYPHI } from '../genotyphi/main';
include { KLEBORATE } from '../kleborate/main';
include { LEGSTA } from '../legsta/main';
include { LISSERO } from '../lissero/main';
include { MENINGOTYPE } from '../meningotype/main';
include { NGMASTER } from '../ngmaster/main';
include { PASTY } from '../pasty/main';
include { PBPTYPER } from '../pbptyper/main';
include { SEQSERO2 } from '../seqsero2/main';
include { SEROBA } from '../seroba/main';
include { SHIGATYPER } from '../shigatyper/main';
include { SHIGEIFINDER } from '../shigeifinder/main';
include { SISTR } from '../sistr/main';
include { SSUISSERO } from '../ssuissero/main';
include { STAPHTYPER } from '../staphtyper/main';
include { TBPROFILER } from '../tbprofiler/main';

workflow MERLIN {
    take:
    assembly // channel: [ val(meta), [ assembly ] ]

    main:
    ch_versions = Channel.empty()

    // ID potential species
    MERLINDIST(assembly)

    // Escherichia/Shigella
    MERLINDIST.out.escherichia.map{meta, assembly, found -> [meta, assembly]}.set{ ch_escherichia }
    MERLINDIST.out.escherichia_fq.map{meta, reads, found -> [meta, reads]}.set{ ch_escherichia_fq }
    ECTYPER(ch_escherichia)
    ch_versions = ch_versions.mix(ECTYPER.out.versions.first())
    SHIGATYPER(ch_escherichia_fq)
    ch_versions = ch_versions.mix(SHIGATYPER.out.versions.first())
    SHIGEIFINDER(ch_escherichia)
    ch_versions = ch_versions.mix(SHIGEIFINDER.out.versions.first())

    // Haemophilus
    MERLINDIST.out.haemophilus.map{meta, assembly, found -> [meta, assembly]}.set{ ch_haemophilus }
    HICAP(ch_haemophilus)
    ch_versions = ch_versions.mix(HICAP.out.versions.first())
    HPSUISSERO(ch_haemophilus)
    ch_versions = ch_versions.mix(HPSUISSERO.out.versions.first())

    // Klebsiella
    MERLINDIST.out.klebsiella.map{meta, assembly, found -> [meta, assembly]}.set{ ch_klebsiella }
    KLEBORATE(ch_klebsiella)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())

    // Legionella 
    MERLINDIST.out.legionella.map{meta, assembly, found -> [meta, assembly]}.set{ ch_legionella }
    LEGSTA(ch_legionella)
    ch_versions = ch_versions.mix(LEGSTA.out.versions.first())

    // Listeria 
    MERLINDIST.out.listeria.map{meta, assembly, found -> [meta, assembly]}.set{ ch_listeria }
    LISSERO(ch_listeria)
    ch_versions = ch_versions.mix(LISSERO.out.versions.first())

    // Mycobacterium 
    MERLINDIST.out.mycobacterium_fq.map{meta, reads, found -> [meta, reads]}.set{ ch_mycobacterium_fq }
    TBPROFILER(ch_mycobacterium_fq)
    ch_versions = ch_versions.mix(TBPROFILER.out.versions.first())

    // Neisseria 
    MERLINDIST.out.neisseria.map{meta, assembly, found -> [meta, assembly]}.set{ ch_neisseria }
    MENINGOTYPE(ch_neisseria)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions.first())
    NGMASTER(ch_neisseria)
    ch_versions = ch_versions.mix(NGMASTER.out.versions.first())

    // Pseudomonas 
    MERLINDIST.out.pseudomonas.map{meta, assembly, found -> [meta, assembly]}.set{ ch_pseudomonas }
    PASTY(ch_pseudomonas)
    ch_versions = ch_versions.mix(PASTY.out.versions.first())

    // Salmonella 
    MERLINDIST.out.salmonella.map{meta, assembly, found -> [meta, assembly]}.set{ ch_salmonella }
    MERLINDIST.out.salmonella_fq.map{meta, reads, found -> [meta, reads]}.set{ ch_salmonella_fq }
    GENOTYPHI(ch_salmonella_fq)
    ch_versions = ch_versions.mix(GENOTYPHI.out.versions.first())
    SEQSERO2(ch_salmonella)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions.first())
    SISTR(ch_salmonella)
    ch_versions = ch_versions.mix(SISTR.out.versions.first())

    // Staphylococcus 
    MERLINDIST.out.staphylococcus.map{meta, assembly, found -> [meta, assembly]}.set{ ch_staphylococcus }
    STAPHTYPER(ch_staphylococcus)
    ch_versions = ch_versions.mix(STAPHTYPER.out.versions)

    // Streptococcus 
    MERLINDIST.out.streptococcus.map{meta, assembly, found -> [meta, assembly]}.set{ ch_streptococcus }
    MERLINDIST.out.streptococcus_fq.map{meta, reads, found -> [meta, reads]}.set{ ch_streptococcus_fq }
    EMMTYPER(ch_streptococcus)
    ch_versions = ch_versions.mix(EMMTYPER.out.versions.first())
    PBPTYPER(ch_streptococcus)
    ch_versions = ch_versions.mix(PBPTYPER.out.versions.first())
    SEROBA(ch_streptococcus_fq)
    ch_versions = ch_versions.mix(SEROBA.out.versions.first())
    SSUISSERO(ch_streptococcus)
    ch_versions = ch_versions.mix(SSUISSERO.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
