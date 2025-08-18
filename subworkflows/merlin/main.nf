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

    // ID potential species
    MERLINDIST(assembly, mash_db)
    ch_versions = ch_versions.mix(MERLINDIST.out.versions)
    ch_logs = ch_logs.mix(MERLINDIST.out.logs)

    // Escherichia/Shigella
    MERLINDIST.out.escherichia.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_escherichia }
    MERLINDIST.out.escherichia_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_escherichia_fq }
    MERLINDIST.out.escherichia_fna_fq.map{_meta, _assembly, _reads, _found -> [_meta, _assembly, _reads]}.set{ ch_escherichia_fna_fq }
    CLERMONTYPING(ch_escherichia)
    ch_versions = ch_versions.mix(CLERMONTYPING.out.versions.first())
    ch_logs = ch_logs.mix(CLERMONTYPING.out.logs)

    ECTYPER(ch_escherichia)
    ch_versions = ch_versions.mix(ECTYPER.out.versions.first())
    ch_logs = ch_logs.mix(ECTYPER.out.logs)

    SHIGAPASS(ch_escherichia)
    ch_versions = ch_versions.mix(SHIGAPASS.out.versions.first())
    ch_logs = ch_logs.mix(SHIGAPASS.out.logs)

    SHIGATYPER(ch_escherichia_fq)
    ch_versions = ch_versions.mix(SHIGATYPER.out.versions.first())
    ch_logs = ch_logs.mix(SHIGATYPER.out.logs)

    SHIGEIFINDER(ch_escherichia)
    ch_versions = ch_versions.mix(SHIGEIFINDER.out.versions.first())
    ch_logs = ch_logs.mix(SHIGEIFINDER.out.logs)

    STECFINDER(ch_escherichia_fna_fq)
    ch_versions = ch_versions.mix(STECFINDER.out.versions.first())
    ch_logs = ch_logs.mix(STECFINDER.out.logs)

    // Haemophilus
    MERLINDIST.out.haemophilus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_haemophilus }
    HICAP(ch_haemophilus)
    ch_versions = ch_versions.mix(HICAP.out.versions.first())
    ch_logs = ch_logs.mix(HICAP.out.logs)

    HPSUISSERO(ch_haemophilus)
    ch_versions = ch_versions.mix(HPSUISSERO.out.versions.first())
    ch_logs = ch_logs.mix(HPSUISSERO.out.logs)

    // Klebsiella
    MERLINDIST.out.klebsiella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_klebsiella }
    KLEBORATE(ch_klebsiella)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())
    ch_logs = ch_logs.mix(KLEBORATE.out.logs)

    // Legionella 
    MERLINDIST.out.legionella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_legionella }
    LEGSTA(ch_legionella)
    ch_versions = ch_versions.mix(LEGSTA.out.versions.first())
    ch_logs = ch_logs.mix(LEGSTA.out.logs)

    // Listeria 
    MERLINDIST.out.listeria.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_listeria }
    LISSERO(ch_listeria)
    ch_versions = ch_versions.mix(LISSERO.out.versions.first())
    ch_logs = ch_logs.mix(LISSERO.out.logs)

    // Mycobacterium 
    MERLINDIST.out.mycobacterium_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_mycobacterium_fq }
    TBPROFILER(ch_mycobacterium_fq)
    ch_versions = ch_versions.mix(TBPROFILER.out.versions.first())
    ch_logs = ch_logs.mix(TBPROFILER.out.logs)

    // Neisseria 
    MERLINDIST.out.neisseria.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_neisseria }
    MENINGOTYPE(ch_neisseria)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions.first())
    ch_logs = ch_logs.mix(MENINGOTYPE.out.logs)

    NGMASTER(ch_neisseria)
    ch_versions = ch_versions.mix(NGMASTER.out.versions.first())
    ch_logs = ch_logs.mix(NGMASTER.out.logs)

    // Pseudomonas 
    MERLINDIST.out.pseudomonas.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_pseudomonas }
    PASTY(ch_pseudomonas)
    ch_versions = ch_versions.mix(PASTY.out.versions.first())
    ch_logs = ch_logs.mix(PASTY.out.logs)

    // Salmonella 
    MERLINDIST.out.salmonella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_salmonella }
    MERLINDIST.out.salmonella_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_salmonella_fq }
    GENOTYPHI(ch_salmonella_fq)
    ch_versions = ch_versions.mix(GENOTYPHI.out.versions.first())
    ch_logs = ch_logs.mix(GENOTYPHI.out.logs)

    SEQSERO2(ch_salmonella)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions.first())
    ch_logs = ch_logs.mix(SEQSERO2.out.logs)

    SISTR(ch_salmonella)
    ch_versions = ch_versions.mix(SISTR.out.versions.first())
    ch_logs = ch_logs.mix(SISTR.out.logs)

    // Staphylococcus 
    MERLINDIST.out.staphylococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_staphylococcus }
    STAPHTYPER(ch_staphylococcus)
    ch_versions = ch_versions.mix(STAPHTYPER.out.versions)
    ch_logs = ch_logs.mix(STAPHTYPER.out.logs)

    // Streptococcus 
    MERLINDIST.out.streptococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_streptococcus }
    MERLINDIST.out.streptococcus_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_streptococcus_fq }
    EMMTYPER(ch_streptococcus)
    ch_versions = ch_versions.mix(EMMTYPER.out.versions.first())
    ch_logs = ch_logs.mix(EMMTYPER.out.logs)

    PBPTYPER(ch_streptococcus)
    ch_versions = ch_versions.mix(PBPTYPER.out.versions.first())
    ch_logs = ch_logs.mix(PBPTYPER.out.logs)

    SEROBA(ch_streptococcus_fq)
    ch_versions = ch_versions.mix(SEROBA.out.versions.first())
    ch_logs = ch_logs.mix(SEROBA.out.logs)

    SSUISSERO(ch_streptococcus)
    ch_versions = ch_versions.mix(SSUISSERO.out.versions.first())
    ch_logs = ch_logs.mix(SSUISSERO.out.logs)

    emit:
    logs = ch_logs
    nf_logs = CLERMONTYPING.out.nf_logs.mix(
        ECTYPER.out.nf_logs,
        EMMTYPER.out.nf_logs,
        GENOTYPHI.out.nf_logs,
        HICAP.out.nf_logs,
        HPSUISSERO.out.nf_logs,
        KLEBORATE.out.nf_logs,
        LEGSTA.out.nf_logs,
        LISSERO.out.nf_logs,
        MENINGOTYPE.out.nf_logs,
        MERLINDIST.out.nf_logs,
        NGMASTER.out.nf_logs,
        PASTY.out.nf_logs,
        PBPTYPER.out.nf_logs,
        SEQSERO2.out.nf_logs,
        SEROBA.out.nf_logs,
        SHIGAPASS.out.nf_logs,
        SHIGATYPER.out.nf_logs,
        SHIGEIFINDER.out.nf_logs,
    )
    versions = ch_versions
}
