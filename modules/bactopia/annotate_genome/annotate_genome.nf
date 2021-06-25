nextflow.enable.dsl = 2

process ANNOTATE_GENOME {
    /* Annotate the assembly using Prokka, use a proteins FASTA if available */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "annotation/${sample}*"

    input:
    tuple val(sample), val(single_end), file(fq), file(fasta), file(total_contigs)
    file prokka_proteins
    file prodigal_tf

    output:
    file "annotation/${sample}*"
    tuple val(sample), file("annotation/${sample}.{ffn,ffn.gz}"),emit: PLASMID_BLAST,optional: true
    tuple val(sample),
        file("annotation/${sample}.{ffn,ffn.gz}"),
        file("annotation/${sample}.{faa,faa.gz}"),emit: ANTIMICROBIAL_RESISTANCE, optional: true
    file "${task.process}/*" optional true

    shell:
    gunzip_fasta = fasta.getName().replace('.gz', '')
    contig_count = total_contigs.getName().replace('total_contigs_', '')
    genus = "Genus"
    species = "species"
    proteins = ""
    if (prokka_proteins.getName() != 'EMPTY_PROTEINS') {
        proteins = "--proteins ${prokka_proteins}"
        if (SPECIES.contains("-")) {
            genus = SPECIES.split('-')[0].capitalize()
            species = SPECIES.split('-')[1]
        } else {
            genus = SPECIES.capitalize()
            species = "spp."
        }
    }

    prodigal = ""
    if (prodigal_tf.getName() != 'EMPTY_TF' && !params.skip_prodigal_tf) {
        prodigal = "--prodigaltf ${prodigal_tf}"
    }

    compliant = params.compliant ? "--compliant" : ""
    locustag = "--locustag ${sample}"
    renamed = false
    // Contig ID must <= 37 characters
    if ("gnl|${params.centre}|${sample}_${contig_count}".length() > 37) {
        locustag = ""
        compliant = "--compliant"
        renamed = true
    }
    addgenes = params.nogenes ? "" : "--addgenes"
    addmrna = params.addmrna ? "--addmrna" : ""
    rawproduct = params.rawproduct ? "--rawproduct" : ""
    cdsrnaolap = params.cdsrnaolap ? "--cdsrnaolap" : ""
    norrna = params.norrna ? "--norrna" : ""
    notrna = params.notrna ? "--notrna" : ""
    rnammer = params.rnammer ? "--rnammer" : ""
    rfam = params.rnammer ? "--rfam" : ""
    '''
    LOG_DIR="!{task.process}"
    mkdir -p ${LOG_DIR}/

    # Print captured STDERR incase of exit
    function print_stderr {
        cat .command.err 1>&2
        ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
    }
    trap print_stderr EXIT

    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
    if [[ !{params.compress} == "true" ]]; then
        gunzip -f !{fasta}
    fi

    if [ "!{renamed}" == "true" ]; then
        echo "Original sample name (!{sample}) not used due to creating a contig ID >37 characters"
    fi

    # Verify AWS files were staged
    if [[ ! -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "true" ]; then
            check-staging.py --fq1 !{fq[0]} --assembly !{gunzip_fasta} --is_single
        else
            check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --assembly !{gunzip_fasta}
        fi
    fi

    # Prokka Version
    echo "# Prokka Version" >> ${LOG_DIR}/!{task.process}.versions
    prokka --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
    prokka --outdir annotation \
        --force \
        --prefix '!{sample}' \
        --genus '!{genus}' \
        --species '!{species}' \
        --evalue '!{params.prokka_evalue}' \
        --coverage !{params.prokka_coverage} \
        --cpus !{task.cpus} \
        --centre '!{params.centre}' \
        --mincontiglen !{params.min_contig_len} \
        !{locustag} \
        !{prodigal} \
        !{addgenes} \
        !{compliant} \
        !{proteins} \
        !{rawproduct} \
        !{cdsrnaolap} \
        !{addmrna} \
        !{norrna} \
        !{notrna} \
        !{rnammer} \
        !{rfam} \
        !{gunzip_fasta} > ${LOG_DIR}/prokka.out 2> ${LOG_DIR}/prokka.err

    if [[ !{params.compress} == "true" ]]; then
        find annotation/ -type f -not -name "*.txt" -and -not -name "*.log*" | \
            xargs -I {} pigz -n --best -p !{task.cpus} {}
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
        cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir annotation
    mkdir ${task.process}
    touch annotation/${sample}
    touch annotation/${sample}.ffn
    touch annotation/${sample}.ffn.gz
    touch annotation/${sample}.faa
    touch annotation/${sample}.faa.gz
    touch "${task.process}/${sample}"
    """
}


//###############
//Module testing
//###############

workflow test{
    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        file(params.fq),
        file(params.fasta),
        file(params.total_contigs)
        ])
    TEST_PARAMS_CH2 = Channel.of(
        file(params.prokka_proteins)
        )
    TEST_PARAMS_CH3 = Channel.of(
        file(params.prodigal_tf)
        )

    annotate_genome(TEST_PARAMS_CH,TEST_PARAMS_CH2,TEST_PARAMS_CH3)
}
