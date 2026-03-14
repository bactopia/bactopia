/**
 * MinER assisted species-specific bactopia tool seLectIoN.
 *
 * This subworkflow performs intelligent species identification and selects appropriate
 * species-specific typing tools based on the detected organism. It first identifies
 * potential species using MinHash distance estimation, then runs species-specific
 * subworkflows for detailed characterization including serotyping, MLST, virulence
 * factor detection, and antimicrobial resistance profiling.
 *
 * @status stable
 * @keywords species, identification, typing, serotype, virulence
 * @tags complexity:complex input-type:single output-type:multiple features:conditional-logic, components
 * @citation mash
 *
 * @subworkflows merlindist, clermontyping, ectyper, emmtyper, genotyphi, hicap, hpsuissero, kleborate,
 *              legsta, lissero, meningotype, ngmaster, pasty, pbptyper, seqsero2, seroba, shigapass,
 *              shigatyper, shigeifinder, sistr, ssuissero, staphtyper, stecfinder, tbprofiler
 *
 * @input record(meta, fna, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembly file for species identification and typing
 * - `r1`: Illumina R1 reads (paired-end) or null
 * - `r2`: Illumina R2 reads (paired-end) or null
 * - `se`: Single-end Illumina reads or null
 * - `lr`: Long reads (ONT/PacBio) or null
 *
 * @input mash_db
 * Mash sketch database for rapid species identification
 *
 * @input emmtyper_blastdb
 * EMMTyper BLAST database for Streptococcus pyogenes emm typing (optional)
 *
 * @input hicap_database_dir
 * HiCAP database directory for Haemophilus influenzae serotyping (optional)
 *
 * @input hicap_model_fp
 * HiCAP HMM model file for improved detection (optional)
 *
 * @input staphtyper_repeats
 * Staphylococcus aureus repeat sequences for spa typing (optional)
 *
 * @input staphtyper_repeat_order
 * Staphylococcus aureus repeat order file for spa typing (optional)
 *
 * @output dist_outputs
 * - `dist`: The raw Mash distance results
 * - `fna`: Passthrough of assembled contigs
 * - `r1`: Passthrough of Illumina R1 reads
 * - `r2`: Passthrough of Illumina R2 reads
 * - `se`: Passthrough of single-end reads
 * - `lr`: Passthrough of long reads
 * - `escherichia`: Conditional marker file triggering Escherichia analysis tools
 * - `haemophilus`: Conditional marker file triggering Haemophilus analysis tools
 * - `klebsiella`: Conditional marker file triggering Klebsiella analysis tools
 * - `legionella`: Conditional marker file triggering Legionella analysis tools
 * - `listeria`: Conditional marker file triggering Listeria analysis tools
 * - `mycobacterium`: Conditional marker file triggering Mycobacterium analysis tools
 * - `neisseria`: Conditional marker file triggering Neisseria analysis tools
 * - `pseudomonas`: Conditional marker file triggering Pseudomonas analysis tools
 * - `salmonella`: Conditional marker file triggering Salmonella analysis tools
 * - `staphylococcus`: Conditional marker file triggering Staphylococcus analysis tools
 * - `streptococcus`: Conditional marker file triggering Streptococcus analysis tools
 * - `genus`: A marker file indicating the detected genus (for debugging)
 *
 * @output sample_outputs
 */
nextflow.preview.types = true

include { MERLINDIST    } from '../merlindist/main';
include { CLERMONTYPING } from '../clermontyping/main';
include { ECTYPER       } from '../ectyper/main';
include { EMMTYPER      } from '../emmtyper/main';
include { HICAP         } from '../hicap/main';
include { HPSUISSERO    } from '../hpsuissero/main';
include { GENOTYPHI     } from '../genotyphi/main';
include { KLEBORATE     } from '../kleborate/main';
include { LEGSTA        } from '../legsta/main';
include { LISSERO       } from '../lissero/main';
include { MENINGOTYPE   } from '../meningotype/main';
include { NGMASTER      } from '../ngmaster/main';
include { PASTY         } from '../pasty/main';
include { PBPTYPER      } from '../pbptyper/main';
include { SEQSERO2      } from '../seqsero2/main';
include { SEROBA        } from '../seroba/main';
include { SHIGAPASS     } from '../shigapass/main';
include { SHIGATYPER    } from '../shigatyper/main';
include { SHIGEIFINDER  } from '../shigeifinder/main';
include { SISTR         } from '../sistr/main';
include { SSUISSERO     } from '../ssuissero/main';
include { STAPHTYPER    } from '../staphtyper/main';
include { STECFINDER    } from '../stecfinder/main';
include { TBPROFILER    } from '../tbprofiler/main';

workflow MERLIN {
    take:
    assembly: Channel<Record>
    mash_db: Path
    emmtyper_blastdb: Path?
    hicap_database_dir: Path?
    hicap_model_fp: Path?
    staphtyper_repeats: Path?
    staphtyper_repeat_order: Path?

    main:
    // ID potential species
    MERLINDIST(assembly, mash_db)

    // Escherichia/Shigella
    // Assembly-only: filter by genus marker, wrap fna in list for Set<Path> compatibility
    ch_escherichia = MERLINDIST.out.sample_outputs
        .filter { r -> r.escherichia != null }
        .map { r -> [r.meta, [r.fna]] }
    // Reads-only: collect non-null reads into list for Set<Path> compatibility
    ch_escherichia_fq = MERLINDIST.out.sample_outputs
        .filter { r -> r.escherichia != null }
        .map { r ->
            def reads = [r.r1, r.r2, r.se, r.lr].findAll { it != null }
            [r.meta, reads]
        }
    // Assembly + reads: pass as 6-slot tuple for STECFINDER
    ch_escherichia_fna_fq = MERLINDIST.out.sample_outputs
        .filter { r -> r.escherichia != null }
        .map { r -> [r.meta, r.fna, r.r1, r.r2, r.se, r.lr] }
    CLERMONTYPING(ch_escherichia)
    ECTYPER(ch_escherichia)
    SHIGAPASS(ch_escherichia)
    SHIGATYPER(ch_escherichia_fq)
    SHIGEIFINDER(ch_escherichia)
    STECFINDER(ch_escherichia_fna_fq)

    // Haemophilus
    ch_haemophilus = MERLINDIST.out.sample_outputs
        .filter { r -> r.haemophilus != null }
        .map { r -> [r.meta, [r.fna]] }
    HICAP(ch_haemophilus, hicap_database_dir, hicap_model_fp)
    HPSUISSERO(ch_haemophilus)

    // Klebsiella
    ch_klebsiella = MERLINDIST.out.sample_outputs
        .filter { r -> r.klebsiella != null }
        .map { r -> [r.meta, [r.fna]] }
    KLEBORATE(ch_klebsiella)

    // Legionella
    ch_legionella = MERLINDIST.out.sample_outputs
        .filter { r -> r.legionella != null }
        .map { r -> [r.meta, [r.fna]] }
    LEGSTA(ch_legionella)

    // Listeria
    ch_listeria = MERLINDIST.out.sample_outputs
        .filter { r -> r.listeria != null }
        .map { r -> [r.meta, [r.fna]] }
    LISSERO(ch_listeria)

    // Mycobacterium
    ch_mycobacterium_fq = MERLINDIST.out.sample_outputs
        .filter { r -> r.mycobacterium != null }
        .map { r ->
            def reads = [r.r1, r.r2, r.se, r.lr].findAll { it != null }
            [r.meta, reads]
        }
    TBPROFILER(ch_mycobacterium_fq)

    // Neisseria
    ch_neisseria = MERLINDIST.out.sample_outputs
        .filter { r -> r.neisseria != null }
        .map { r -> [r.meta, [r.fna]] }
    MENINGOTYPE(ch_neisseria)
    NGMASTER(ch_neisseria)

    // Pseudomonas
    ch_pseudomonas = MERLINDIST.out.sample_outputs
        .filter { r -> r.pseudomonas != null }
        .map { r -> [r.meta, [r.fna]] }
    PASTY(ch_pseudomonas)

    // Salmonella
    ch_salmonella = MERLINDIST.out.sample_outputs
        .filter { r -> r.salmonella != null }
        .map { r -> [r.meta, [r.fna]] }
    ch_salmonella_fq = MERLINDIST.out.sample_outputs
        .filter { r -> r.salmonella != null }
        .map { r ->
            def reads = [r.r1, r.r2, r.se, r.lr].findAll { it != null }
            [r.meta, reads]
        }
    GENOTYPHI(ch_salmonella_fq)
    SEQSERO2(ch_salmonella)
    SISTR(ch_salmonella)

    // Staphylococcus
    ch_staphylococcus = MERLINDIST.out.sample_outputs
        .filter { r -> r.staphylococcus != null }
        .map { r -> [r.meta, [r.fna]] }
    STAPHTYPER(ch_staphylococcus, staphtyper_repeats, staphtyper_repeat_order)

    // Streptococcus
    ch_streptococcus = MERLINDIST.out.sample_outputs
        .filter { r -> r.streptococcus != null }
        .map { r -> [r.meta, [r.fna]] }
    ch_streptococcus_fq = MERLINDIST.out.sample_outputs
        .filter { r -> r.streptococcus != null }
        .map { r ->
            def reads = [r.r1, r.r2, r.se, r.lr].findAll { it != null }
            [r.meta, reads]
        }
    EMMTYPER(ch_streptococcus, emmtyper_blastdb)
    PBPTYPER(ch_streptococcus)
    SEROBA(ch_streptococcus_fq)
    SSUISSERO(ch_streptococcus)

    emit:
    dist_outputs = MERLINDIST.out.sample_outputs
    sample_outputs = CLERMONTYPING.out.sample_outputs.mix(
        ECTYPER.out.sample_outputs,
        EMMTYPER.out.sample_outputs,
        GENOTYPHI.out.sample_outputs,
        HICAP.out.sample_outputs,
        HPSUISSERO.out.sample_outputs,
        KLEBORATE.out.sample_outputs,
        LEGSTA.out.sample_outputs,
        LISSERO.out.sample_outputs,
        MENINGOTYPE.out.sample_outputs,
        NGMASTER.out.sample_outputs,
        PASTY.out.sample_outputs,
        PBPTYPER.out.sample_outputs,
        SEQSERO2.out.sample_outputs,
        SEROBA.out.sample_outputs,
        SHIGAPASS.out.sample_outputs,
        SHIGATYPER.out.sample_outputs,
        SHIGEIFINDER.out.sample_outputs,
        STECFINDER.out.sample_outputs,
        SISTR.out.sample_outputs,
        SSUISSERO.out.sample_outputs,
        STAPHTYPER.out.sample_outputs,
        TBPROFILER.out.sample_outputs
    )
}
