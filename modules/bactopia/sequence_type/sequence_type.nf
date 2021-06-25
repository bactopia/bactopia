nextflow.enable.dsl = 2

process SEQUENCE_TYPE {
    /* Determine MLST types using ARIBA and BLAST */
    tag "${sample} - ${schema} - ${method}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/mlst/${schema}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${method}/*"

    input:
    tuple val(sample), val(single_end), path(fq), path(assembly)
    each path(dataset)

    output:
    file "${method}/*"
    file "${task.process}/*" optional true

    when:
    MLST_DATABASES.isEmpty() == false

    shell:
    method = dataset =~ /.*blastdb.*/ ? 'blast' : 'ariba'
    dataset_tarball = path(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '').split('-')[1]
    schema = dataset_tarball.split('-')[0]
    noclean = params.ariba_no_clean ? "--noclean" : ""
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    '''
    LOG_DIR="!{task.process}"
    tar -xzvf !{dataset_tarball}
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}-!{method}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}-!{method}.versions

    # Verify AWS files were staged
    if [[ ! -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "true" ]; then
            check-staging.py --fq1 !{fq[0]} --assembly !{assembly} --is_single
        else
            check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --assembly !{assembly}
        fi
    fi

    if [ "!{method}" == "blast" ]; then
        echo "# mlst-blast.py Version" >> ${LOG_DIR}/!{task.process}-!{method}.versions
        mlst-blast.py --version >> ${LOG_DIR}/!{task.process}-!{method}.versions 2>&1
        mkdir -p blast
        if [[ !{params.compress} == "true" ]]; then
            mlst-blast.py !{assembly} !{dataset_name} blast/!{sample}-blast.json \
                --cpu !{task.cpus} --compressed
        else
            mlst-blast.py !{assembly} !{dataset_name} blast/!{sample}-blast.json \
                --cpu !{task.cpus}
        fi
    elif [ "!{method}" == "ariba" ]; then
        if [ "!{single_end}" == "false" ]; then
            echo "# Ariba Version" >> ${LOG_DIR}/!{task.process}-!{method}.versions
            ariba version >> ${LOG_DIR}/!{task.process}-!{method}.versions 2>&1
            mv !{dataset_name}/ref_db ./
            ariba run ref_db !{fq[0]} !{fq[1]} ariba \
                --nucmer_min_id !{params.nucmer_min_id} \
                --nucmer_min_len !{params.nucmer_min_len} \
                --nucmer_breaklen !{params.nucmer_breaklen} \
                --assembly_cov !{params.assembly_cov} \
                --min_scaff_depth !{params.min_scaff_depth} \
                --assembled_threshold !{params.assembled_threshold} \
                --gene_nt_extend !{params.gene_nt_extend} \
                --unique_threshold !{params.unique_threshold} \
                --threads !{task.cpus} \
                --force \
                --verbose !{noclean} !{spades_options}
        else
            mkdir -p ariba
            echo "Ariba cannot be run on single end reads" > ariba/ariba-not-run.txt
        fi
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}-!{method}.err
        cp .command.out ${LOG_DIR}/!{task.process}-!{method}.out
        cp .command.sh ${LOG_DIR}/!{task.process}-!{method}.sh  || :
        cp .command.trace ${LOG_DIR}/!{task.process}-!{method}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    method = dataset =~ /.*blastdb.*/ ? 'blast' : 'ariba'
    dataset_tarball = path(dataset).getName()
    schema = dataset_tarball.split('-')[0]
    """
    mkdir ${method}
    mkdir ${task.process}
    touch ${method}/${sample}
    touch ${task.process}/${sample}
    """
}

//###############
//Module testing
//###############

workflow test{

    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        path(params.fq),
        path(params.assembly)
        ])
    TEST_PARAMS_CH2 = Channel.of(
        path(params.dataset_blast)
        path(params.dataset_ariba))

    sequence_type(TEST_PARAMS_CH,TEST_PARAMS_CH2.collect())
}
