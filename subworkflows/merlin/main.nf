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
    // ID potential species
    MERLINDIST(assembly, mash_db)

    // Escherichia/Shigella
    MERLINDIST.out.escherichia.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_escherichia }
    MERLINDIST.out.escherichia_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_escherichia_fq }
    MERLINDIST.out.escherichia_fna_fq.map{_meta, _assembly, _reads, _found -> [_meta, _assembly, _reads]}.set{ ch_escherichia_fna_fq }
    CLERMONTYPING(ch_escherichia)
    ECTYPER(ch_escherichia)
    SHIGAPASS(ch_escherichia)
    SHIGATYPER(ch_escherichia_fq)
    SHIGEIFINDER(ch_escherichia)
    STECFINDER(ch_escherichia_fna_fq)

    // Haemophilus
    MERLINDIST.out.haemophilus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_haemophilus }
    HICAP(ch_haemophilus)
    HPSUISSERO(ch_haemophilus)

    // Klebsiella
    MERLINDIST.out.klebsiella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_klebsiella }
    KLEBORATE(ch_klebsiella)

    // Legionella 
    MERLINDIST.out.legionella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_legionella }
    LEGSTA(ch_legionella)

    // Listeria 
    MERLINDIST.out.listeria.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_listeria }
    LISSERO(ch_listeria)

    // Mycobacterium 
    MERLINDIST.out.mycobacterium_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_mycobacterium_fq }
    TBPROFILER(ch_mycobacterium_fq)

    // Neisseria 
    MERLINDIST.out.neisseria.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_neisseria }
    MENINGOTYPE(ch_neisseria)
    NGMASTER(ch_neisseria)

    // Pseudomonas 
    MERLINDIST.out.pseudomonas.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_pseudomonas }
    PASTY(ch_pseudomonas)

    // Salmonella 
    MERLINDIST.out.salmonella.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_salmonella }
    MERLINDIST.out.salmonella_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_salmonella_fq }
    GENOTYPHI(ch_salmonella_fq)
    SEQSERO2(ch_salmonella)
    SISTR(ch_salmonella)

    // Staphylococcus 
    MERLINDIST.out.staphylococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_staphylococcus }
    STAPHTYPER(ch_staphylococcus)

    // Streptococcus 
    MERLINDIST.out.streptococcus.map{_meta, _assembly, _found -> [_meta, _assembly]}.set{ ch_streptococcus }
    MERLINDIST.out.streptococcus_fq.map{_meta, _reads, _found -> [_meta, _reads]}.set{ ch_streptococcus_fq }
    EMMTYPER(ch_streptococcus)
    PBPTYPER(ch_streptococcus)
    SEROBA(ch_streptococcus_fq)
    SSUISSERO(ch_streptococcus)

    emit:
    results = MERLINDIST.out.results.mix(
        CLERMONTYPING.out.results,
        ECTYPER.out.results,
        EMMTYPER.out.results,
        GENOTYPHI.out.results,
        HICAP.out.results,
        HPSUISSERO.out.results,
        KLEBORATE.out.results,
        LEGSTA.out.results,
        LISSERO.out.results,
        MENINGOTYPE.out.results,
        NGMASTER.out.results,
        PASTY.out.results,
        PBPTYPER.out.results,
        SEQSERO2.out.results,
        SEROBA.out.results,
        SHIGAPASS.out.results,
        SHIGATYPER.out.results,
        SHIGEIFINDER.out.results,
        STECFINDER.out.results,
        SISTR.out.results,
        SSUISSERO.out.results,
        STAPHTYPER.out.results,
        TBPROFILER.out.results
    )
    logs = MERLINDIST.out.logs.mix(
        CLERMONTYPING.out.logs,
        ECTYPER.out.logs,
        EMMTYPER.out.logs,
        GENOTYPHI.out.logs,
        HICAP.out.logs,
        HPSUISSERO.out.logs,
        KLEBORATE.out.logs,
        LEGSTA.out.logs,
        LISSERO.out.logs,
        MENINGOTYPE.out.logs,
        NGMASTER.out.logs,
        PASTY.out.logs,
        PBPTYPER.out.logs,
        SEQSERO2.out.logs,
        SEROBA.out.logs,
        SHIGAPASS.out.logs,
        SHIGATYPER.out.logs,
        SHIGEIFINDER.out.logs,
        STECFINDER.out.logs,
        SISTR.out.logs,
        SSUISSERO.out.logs,
        STAPHTYPER.out.logs,
        TBPROFILER.out.logs
    )
    nf_logs = MERLINDIST.out.nf_logs.mix(
        CLERMONTYPING.out.nf_logs,
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
        STECFINDER.out.nf_logs,
        SISTR.out.nf_logs,
        SSUISSERO.out.nf_logs,
        STAPHTYPER.out.nf_logs,
        TBPROFILER.out.nf_logs
    )
    versions = MERLINDIST.out.versions.mix(
        CLERMONTYPING.out.versions,
        ECTYPER.out.versions,
        EMMTYPER.out.versions,
        GENOTYPHI.out.versions,
        HICAP.out.versions,
        HPSUISSERO.out.versions,
        KLEBORATE.out.versions,
        LEGSTA.out.versions,
        LISSERO.out.versions,
        MENINGOTYPE.out.versions,
        NGMASTER.out.versions,
        PASTY.out.versions,
        PBPTYPER.out.versions,
        SEQSERO2.out.versions,
        SEROBA.out.versions,
        SHIGAPASS.out.versions,
        SHIGATYPER.out.versions,
        SHIGEIFINDER.out.versions,
        STECFINDER.out.versions,
        SISTR.out.versions,
        SSUISSERO.out.versions,
        STAPHTYPER.out.versions,
        TBPROFILER.out.versions
    )
}
