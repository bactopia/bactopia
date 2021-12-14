nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'blast')

process BLAST {
    /*
    Query gene FASTA files against annotated assembly using BLAST
    */
    tag "${meta.id} - ${query}"
    label "blast"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options, logs_subdir: query) }

    input:
    tuple val(meta), path(blastdb, stageAs: 'blastdb/*')
    each path(query)

    output:
    path "results/*"
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    shell:
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
            blastn -db blastdb/!{meta.id} \
                -outfmt 15 \
                -evalue 1 \
                -perc_identity !{params.perc_identity} \
                -qcov_hsp_perc !{params.qcov_hsp_perc} \
                -query - \
                -out temp_json/${name}_{#}.json > blastn.stdout.txt 2> blastn.stderr.txt
        elif [ "!{query}" == "primers" ]; then
            cat ${fasta} | sed -e 's/<[^>]*>//g' |
            parallel --gnu --plain -j !{task.cpus} --recstart '>' -N 1 --pipe \
            blastn -db blastdb/!{meta.id} \
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
            tblastn -db blastdb/!{meta.id} \
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

    cat <<-END_VERSIONS > versions.yml
    blast:
        blastn: $(echo $(blastn -version 2>&1) | sed 's/^.*blastn: //;s/ .*$//')
        parallel: $(echo $(parallel --version 2>&1) | sed 's/^GNU parallel //;s/ .*$//')
        pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
        tblastn: $(echo $(tblastn -version 2>&1) | sed 's/^.*tblastn: //;s/ .*$//') 
    END_VERSIONS
    '''
}
