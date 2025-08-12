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

    emit:
    logs = ch_logs
    nf_logs = CLERMONTYPING.out.nf_begin.mix(
        CLERMONTYPING.out.nf_err,
        CLERMONTYPING.out.nf_log,
        CLERMONTYPING.out.nf_out,
        CLERMONTYPING.out.nf_run,
        CLERMONTYPING.out.nf_sh,
        CLERMONTYPING.out.nf_trace,
        ECTYPER.out.nf_begin,
        ECTYPER.out.nf_err,
        ECTYPER.out.nf_log,
        ECTYPER.out.nf_out,
        ECTYPER.out.nf_run,
        ECTYPER.out.nf_sh,
        ECTYPER.out.nf_trace,
        EMMTYPER.out.nf_begin,
        EMMTYPER.out.nf_err,
        EMMTYPER.out.nf_log,
        EMMTYPER.out.nf_out,
        EMMTYPER.out.nf_run,
        EMMTYPER.out.nf_sh,
        EMMTYPER.out.nf_trace,
        GENOTYPHI.out.nf_begin,
        GENOTYPHI.out.nf_err,
        GENOTYPHI.out.nf_log,
        GENOTYPHI.out.nf_out,
        GENOTYPHI.out.nf_run,
        GENOTYPHI.out.nf_sh,
        GENOTYPHI.out.nf_trace,
        HICAP.out.nf_begin,
        HICAP.out.nf_err,
        HICAP.out.nf_log,
        HICAP.out.nf_out,
        HICAP.out.nf_run,
        HICAP.out.nf_sh,
        HICAP.out.nf_trace,
        HPSUISSERO.out.nf_begin,
        HPSUISSERO.out.nf_err,
        HPSUISSERO.out.nf_log,
        HPSUISSERO.out.nf_out,
        HPSUISSERO.out.nf_run,
        HPSUISSERO.out.nf_sh,
        HPSUISSERO.out.nf_trace,
        KLEBORATE.out.nf_begin,
        KLEBORATE.out.nf_err,
        KLEBORATE.out.nf_log,
        KLEBORATE.out.nf_out,
        KLEBORATE.out.nf_run,
        KLEBORATE.out.nf_sh,
        KLEBORATE.out.nf_trace,
        LEGSTA.out.nf_begin,
        LEGSTA.out.nf_err,
        LEGSTA.out.nf_log,
        LEGSTA.out.nf_out,
        LEGSTA.out.nf_run,
        LEGSTA.out.nf_sh,
        LEGSTA.out.nf_trace,
        LISSERO.out.nf_begin,
        LISSERO.out.nf_err,
        LISSERO.out.nf_log,
        LISSERO.out.nf_out,
        LISSERO.out.nf_run,
        LISSERO.out.nf_sh,
        LISSERO.out.nf_trace,
        MENINGOTYPE.out.nf_begin,
        MENINGOTYPE.out.nf_err,
        MENINGOTYPE.out.nf_log,
        MENINGOTYPE.out.nf_out,
        MENINGOTYPE.out.nf_run,
        MENINGOTYPE.out.nf_sh,
        MENINGOTYPE.out.nf_trace,
        MERLINDIST.out.nf_begin,
        MERLINDIST.out.nf_err,
        MERLINDIST.out.nf_log,
        MERLINDIST.out.nf_out,
        MERLINDIST.out.nf_run,
        MERLINDIST.out.nf_sh,
        MERLINDIST.out.nf_trace,
        NGMASTER.out.nf_begin,
        NGMASTER.out.nf_err,
        NGMASTER.out.nf_log,
        NGMASTER.out.nf_out,
        NGMASTER.out.nf_run,
        NGMASTER.out.nf_sh,
        NGMASTER.out.nf_trace,
        PASTY.out.nf_begin,
        PASTY.out.nf_err,
        PASTY.out.nf_log,
        PASTY.out.nf_out,
        PASTY.out.nf_run,
        PASTY.out.nf_sh,
        PASTY.out.nf_trace,
        PBPTYPER.out.nf_begin,
        PBPTYPER.out.nf_err,
        PBPTYPER.out.nf_log,
        PBPTYPER.out.nf_out,
        PBPTYPER.out.nf_run,
        PBPTYPER.out.nf_sh,
        PBPTYPER.out.nf_trace,
        SEQSERO2.out.nf_begin,
        SEQSERO2.out.nf_err,
        SEQSERO2.out.nf_log,
        SEQSERO2.out.nf_out,
        SEQSERO2.out.nf_run,
        SEQSERO2.out.nf_sh,
        SEQSERO2.out.nf_trace,
        SEROBA.out.nf_begin,
        SEROBA.out.nf_err,
        SEROBA.out.nf_log,
        SEROBA.out.nf_out,
        SEROBA.out.nf_run,
        SEROBA.out.nf_sh,
        SEROBA.out.nf_trace,
        SHIGAPASS.out.nf_begin,
        SHIGAPASS.out.nf_err,
        SHIGAPASS.out.nf_log,
        SHIGAPASS.out.nf_out,
        SHIGAPASS.out.nf_run,
        SHIGAPASS.out.nf_sh,
        SHIGAPASS.out.nf_trace,
        SHIGATYPER.out.nf_begin,
        SHIGATYPER.out.nf_err,
        SHIGATYPER.out.nf_log,
        SHIGATYPER.out.nf_out,
        SHIGATYPER.out.nf_run,
        SHIGATYPER.out.nf_sh,
        SHIGATYPER.out.nf_trace,
        SHIGEIFINDER.out.nf_begin,
        SHIGEIFINDER.out.nf_err,
        SHIGEIFINDER.out.nf_log,
        SHIGEIFINDER.out.nf_out,
        SHIGEIFINDER.out.nf_run,
        SHIGEIFINDER.out.nf_sh,
        SHIGEIFINDER.out.nf_trace,
        SISTR.out.nf_begin,
        SISTR.out.nf_err,
        SISTR.out.nf_log,
        SISTR.out.nf_out,
        SISTR.out.nf_run,
        SISTR.out.nf_sh,
        SISTR.out.nf_trace,
        SSUISSERO.out.nf_begin,
        SSUISSERO.out.nf_err,
        SSUISSERO.out.nf_log,
        SSUISSERO.out.nf_out,
        SSUISSERO.out.nf_run,
        SSUISSERO.out.nf_sh,
        SSUISSERO.out.nf_trace,
        STAPHTYPER.out.nf_begin,
        STAPHTYPER.out.nf_err,
        STAPHTYPER.out.nf_log,
        STAPHTYPER.out.nf_out,
        STAPHTYPER.out.nf_run,
        STAPHTYPER.out.nf_sh,
        STAPHTYPER.out.nf_trace,
        STECFINDER.out.nf_begin,
        STECFINDER.out.nf_err,
        STECFINDER.out.nf_log,
        STECFINDER.out.nf_out,
        STECFINDER.out.nf_run,
        STECFINDER.out.nf_sh,
        STECFINDER.out.nf_trace,
        TBPROFILER.out.nf_begin,
        TBPROFILER.out.nf_err,
        TBPROFILER.out.nf_log,
        TBPROFILER.out.nf_out,
        TBPROFILER.out.nf_run,
        TBPROFILER.out.nf_sh,
        TBPROFILER.out.nf_trace
    )
    versions = ch_versions
}
