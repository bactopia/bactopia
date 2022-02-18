//
// merlin - MinmER assisted species-specific bactopia tool seLectIoN
//

include { MERLINDIST } from '../mashdist/main';
include { ECTYPER } from '../ectyper/main';
include { EMMTYPER } from '../emmtyper/main';
include { HICAP } from '../hicap/main';
include { KLEBORATE } from '../kleborate/main';
include { LEGSTA } from '../legsta/main';
include { LISSERO } from '../lissero/main';
include { MENINGOTYPE } from '../meningotype/main';
include { NGMASTER } from '../ngmaster/main';
include { SEQSERO2 } from '../seqsero2/main';
include { SISTR } from '../sistr/main';
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
    ECTYPER(ch_escherichia)
    ch_versions = ch_versions.mix(ECTYPER.out.versions.first())

    // Haemophilus
    MERLINDIST.out.haemophilus.map{meta, assembly, found -> [meta, assembly]}.set{ ch_haemophilus }
    HICAP(ch_haemophilus)
    ch_versions = ch_versions.mix(HICAP.out.versions.first())

    // Klebsiella
    MERLINDIST.out.klebsiella.map{meta, assembly, found -> [meta, assembly]}.set{ ch_klebsiella }
    KLEBORATE(ch_klebsiella)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())

    // Legionella 
    MERLINDIST.out.legionella.map{meta, assembly, found -> [meta, assembly]}.set{ ch_legionella }
    LEGSTA(ch_listeria)
    ch_versions = ch_versions.mix(LEGSTA.out.versions.first())

    // Listeria 
    MERLINDIST.out.listeria.map{meta, assembly, found -> [meta, assembly]}.set{ ch_listeria }
    LISSERO(ch_listeria)
    ch_versions = ch_versions.mix(LISSERO.out.versions.first())

    // Mycobacterium 
    MERLINDIST.out.mycobacterium_fq.map{meta, reads, found -> [meta, reads]}.set{ ch_mycobacterium }
    TBPROFILER(ch_mycobacterium)
    ch_versions = ch_versions.mix(TBPROFILER.out.versions.first())

    // Neisseria 
    MERLINDIST.out.neisseria.map{meta, assembly, found -> [meta, assembly]}.set{ ch_neisseria }
    MENINGOTYPE(ch_neisseria)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions.first())
    NGMASTER(ch_neisseria)
    ch_versions = ch_versions.mix(NGMASTER.out.versions.first())

    // Salmonella 
    MERLINDIST.out.salmonella.map{meta, assembly, found -> [meta, assembly]}.set{ ch_salmonella }
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
    EMMTYPER(ch_streptococcus)
    ch_versions = ch_versions.mix(EMMTYPER.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
