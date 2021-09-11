nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources;save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "blast"

process BLAST {
    /*
    Query gene FASTA files against annotated assembly using BLAST
    */
    tag "${sample} - ${blast_exe}"
    label "max_cpus"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME, logs_subdir: query) }

    input:
    tuple val(sample), path(blastdb, stageAs: 'blastdb/*')
    each path(query)

    output:
    path "results/*"
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    shell:
    blast_exe = query == "proteins" ? 'tblastn' : 'blastn'
    '''
    OUTDIR=results/!{query}
    mkdir -p ${OUTDIR}
    for fasta in !{query}/*; do
        type=`readlink -f ${fasta}`
        name=$(basename "${fasta%.*}")
        mkdir -p temp_json
        if [ "!{query}" == "genes" ]; then
            cat ${fasta} | sed -e 's/<[^>]*>//g' |
            parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
            blastn -db blastdb/!{sample} \
                -outfmt 15 \
                -evalue 1 \
                -perc_identity !{params.perc_identity} \
                -qcov_hsp_perc !{params.qcov_hsp_perc} \
                -query - \
                -out temp_json/${name}_{#}.json > blastn.stdout.txt 2> blastn.stderr.txt
        elif [ "!{query}" == "primers" ]; then
            cat ${fasta} | sed -e 's/<[^>]*>//g' |
            parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
            blastn -db blastdb/!{sample} \
                -outfmt 15 \
                -task blastn \
                -dust no \
                -word_size 7 \
                -perc_identity !{params.perc_identity} \
                -evalue 1 \
                -query - \
                -out temp_json/${name}_{#}.json > blastn.stdout.txt 2> blastn.stderr.txt
        else
            cat ${fasta} | sed -e 's/<[^>]*>//g' |
            parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
            tblastn -db blastdb/!{sample} \
                -outfmt 15 \
                -evalue 0.0001 \
                -qcov_hsp_perc !{params.qcov_hsp_perc} \
                -query - \
                -out temp_json/${name}_{#}.json > tblastn.stdout.txt 2> tblastn.stderr.txt
        fi

        merge-blast-json.py temp_json > ${OUTDIR}/${name}.json 2> merge-blast-json.stderr.txt
        rm -rf temp_json

        if [[ !{params.skip_compression} == "false" ]]; then
            pigz -n --best -p !{task.cpus} ${OUTDIR}/${name}.json
        fi
    done

    # Capture version
    !{blast_exe} -version > !{blast_exe}.version.txt 2>&1
    parallel --version > parallel.version.txt 2>&1
    '''

    stub:
    """
    mkdir results/${query}
    touch results/${query}/${sample}.json
    touch results/${query}/${sample}.json.gz
    """
}
