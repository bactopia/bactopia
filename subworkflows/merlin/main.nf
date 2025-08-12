//
// merlin - MinmER assisted species-specific bactopia tool seLectIoN
//

include { MERLINDIST } from '../mashdist/main';
include { CLERMONTYPING } from '../clermontyping/main';
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
include { SHIGAPASS } from '../shigapass/main';
include { SHIGATYPER } from '../shigatyper/main';
include { SHIGEIFINDER } from '../shigeifinder/main';
include { SISTR } from '../sistr/main';
include { SSUISSERO } from '../ssuissero/main';
include { STAPHTYPER } from '../staphtyper/main';
include { STECFINDER } from '../stecfinder/main';
include { TBPROFILER } from '../tbprofiler/main';

workflow MERLIN {
    take:
    assembly // channel: [ val(meta), [ assembly ] ]
    mash_db // channel: [ mash_db ]

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()
    ch_nf_logs = Channel.empty()

    // ID potential species
    MERLINDIST(assembly, mash_db)
    ch_versions = ch_versions.mix(MERLINDIST.out.versions)
    ch_logs = ch_logs.mix(MERLINDIST.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MERLINDIST.out.nf_logs)

    // Escherichia/Shigella
    MERLINDIST.out.escherichia.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_escherichia }
    MERLINDIST.out.escherichia_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_escherichia_fq }
    MERLINDIST.out.escherichia_fna_fq.map{_meta, _assembly, _reads, _found -> [_meta, _assembly, _reads]}.set{ ch_escherichia_fna_fq }
    CLERMONTYPING(ch_escherichia)
    ch_versions = ch_versions.mix(CLERMONTYPING.out.versions.first())
    ch_logs = ch_logs.mix(CLERMONTYPING.out.logs)
    ch_nf_logs = ch_nf_logs.mix(CLERMONTYPING.out.nf_logs)
    ECTYPER(ch_escherichia)
    ch_versions = ch_versions.mix(ECTYPER.out.versions.first())
    ch_logs = ch_logs.mix(ECTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(ECTYPER.out.nf_logs)
    SHIGAPASS(ch_escherichia)
    ch_versions = ch_versions.mix(SHIGAPASS.out.versions.first())
    ch_logs = ch_logs.mix(SHIGAPASS.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SHIGAPASS.out.nf_logs)
    SHIGATYPER(ch_escherichia_fq)
    ch_versions = ch_versions.mix(SHIGATYPER.out.versions.first())
    ch_logs = ch_logs.mix(SHIGATYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SHIGATYPER.out.nf_logs)
    SHIGEIFINDER(ch_escherichia)
    ch_versions = ch_versions.mix(SHIGEIFINDER.out.versions.first())
    ch_logs = ch_logs.mix(SHIGEIFINDER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SHIGEIFINDER.out.nf_logs)
    STECFINDER(ch_escherichia_fna_fq)
    ch_versions = ch_versions.mix(STECFINDER.out.versions.first())
    ch_logs = ch_logs.mix(STECFINDER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(STECFINDER.out.nf_logs)

    // Haemophilus
    MERLINDIST.out.haemophilus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_haemophilus }
    HICAP(ch_haemophilus)
    ch_versions = ch_versions.mix(HICAP.out.versions.first())
    ch_logs = ch_logs.mix(HICAP.out.logs)
    ch_nf_logs = ch_nf_logs.mix(HICAP.out.nf_logs)
    HPSUISSERO(ch_haemophilus)
    ch_versions = ch_versions.mix(HPSUISSERO.out.versions.first())
    ch_logs = ch_logs.mix(HPSUISSERO.out.logs)
    ch_nf_logs = ch_nf_logs.mix(HPSUISSERO.out.nf_logs)

    // Klebsiella
    MERLINDIST.out.klebsiella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_klebsiella }
    KLEBORATE(ch_klebsiella)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())
    ch_logs = ch_logs.mix(KLEBORATE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(KLEBORATE.out.nf_logs)

    // Legionella 
    MERLINDIST.out.legionella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_legionella }
    LEGSTA(ch_legionella)
    ch_versions = ch_versions.mix(LEGSTA.out.versions.first())
    ch_logs = ch_logs.mix(LEGSTA.out.logs)
    ch_nf_logs = ch_nf_logs.mix(LEGSTA.out.nf_logs)

    // Listeria 
    MERLINDIST.out.listeria.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_listeria }
    LISSERO(ch_listeria)
    ch_versions = ch_versions.mix(LISSERO.out.versions.first())
    ch_logs = ch_logs.mix(LISSERO.out.logs)
    ch_nf_logs = ch_nf_logs.mix(LISSERO.out.nf_logs)

    // Mycobacterium 
    MERLINDIST.out.mycobacterium_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_mycobacterium_fq }
    TBPROFILER(ch_mycobacterium_fq)
    ch_versions = ch_versions.mix(TBPROFILER.out.versions.first())
    ch_logs = ch_logs.mix(TBPROFILER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(TBPROFILER.out.nf_logs)

    // Neisseria 
    MERLINDIST.out.neisseria.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_neisseria }
    MENINGOTYPE(ch_neisseria)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions.first())
    ch_logs = ch_logs.mix(MENINGOTYPE.out.logs)
    ch_nf_logs = ch_nf_logs.mix(MENINGOTYPE.out.nf_logs)
    NGMASTER(ch_neisseria)
    ch_versions = ch_versions.mix(NGMASTER.out.versions.first())
    ch_logs = ch_logs.mix(NGMASTER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(NGMASTER.out.nf_logs)

    // Pseudomonas 
    MERLINDIST.out.pseudomonas.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_pseudomonas }
    PASTY(ch_pseudomonas)
    ch_versions = ch_versions.mix(PASTY.out.versions.first())
    ch_logs = ch_logs.mix(PASTY.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PASTY.out.nf_logs)

    // Salmonella 
    MERLINDIST.out.salmonella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_salmonella }
    MERLINDIST.out.salmonella_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_salmonella_fq }
    GENOTYPHI(ch_salmonella_fq)
    ch_versions = ch_versions.mix(GENOTYPHI.out.versions.first())
    ch_logs = ch_logs.mix(GENOTYPHI.out.logs)
    ch_nf_logs = ch_nf_logs.mix(GENOTYPHI.out.nf_logs)
    SEQSERO2(ch_salmonella)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions.first())
    ch_logs = ch_logs.mix(SEQSERO2.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SEQSERO2.out.nf_logs)
    SISTR(ch_salmonella)
    ch_versions = ch_versions.mix(SISTR.out.versions.first())
    ch_logs = ch_logs.mix(SISTR.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SISTR.out.nf_logs)

    // Staphylococcus 
    MERLINDIST.out.staphylococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_staphylococcus }
    STAPHTYPER(ch_staphylococcus)
    ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    ch_logs = ch_logs.mix(STAPHTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(STAPHTYPER.out.nf_logs)

    // Streptococcus 
    MERLINDIST.out.streptococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_streptococcus }
    MERLINDIST.out.streptococcus_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_streptococcus_fq }
    EMMTYPER(ch_streptococcus)
    ch_versions = ch_versions.mix(EMMTYPER.out.versions.first())
    ch_logs = ch_logs.mix(EMMTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(EMMTYPER.out.nf_logs)
    PBPTYPER(ch_streptococcus)
    ch_versions = ch_versions.mix(PBPTYPER.out.versions.first())
    ch_logs = ch_logs.mix(PBPTYPER.out.logs)
    ch_nf_logs = ch_nf_logs.mix(PBPTYPER.out.nf_logs)
    SEROBA(ch_streptococcus_fq)
    ch_versions = ch_versions.mix(SEROBA.out.versions.first())
    ch_logs = ch_logs.mix(SEROBA.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SEROBA.out.nf_logs)
    SSUISSERO(ch_streptococcus)
    ch_versions = ch_versions.mix(SSUISSERO.out.versions.first())
    ch_logs = ch_logs.mix(SSUISSERO.out.logs)
    ch_nf_logs = ch_nf_logs.mix(SSUISSERO.out.nf_logs)

    emit:
    logs = ch_logs
    nf_logs = ch_nf_logs
    versions = ch_versions // channel: [ versions.yml ]
}
