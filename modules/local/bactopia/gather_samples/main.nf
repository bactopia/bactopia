nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options = initOptions(params.containsKey('options') ? params.options : [:], 'gather_samples')
options.is_module = params.wf == 'gather_samples' ? true : false
options.ignore = [".fastq.gz", ".fna.gz"]

process GATHER_SAMPLES {
    /* Gather up input FASTQs for analysis. */
    tag "${meta.id}"
    label "gather_samples"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    input:
    tuple val(meta), path(r1, stageAs: '*???-r1'), path(r2, stageAs: '*???-r2'), path(extra)

    output:
    tuple val(meta), path("fastqs/${meta.id}*.fastq.gz"), path("extra/*.gz"), path("${meta.id}-genome-size.txt"), emit: raw_fastq, optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions
    path "*-{error,merged}.txt", optional: true

    shell:
    meta.original_runtype = meta.runtype
    runtype = meta.original_runtype
    is_assembly = runtype.startsWith('assembly') ? true : false
    is_compressed = extra ? (extra.getName().endsWith('gz') ? true : false) : false
    no_cache = params.no_cache ? '-N' : ''
    archive = params.use_ena ? (task.attempt >= 4 ? "SRA" : "ENA") : "SRA"
    section = runtype == 'assembly_accession' ? (meta.id.startsWith('GCF') ? 'refseq' : 'genbank') : null
    fcov = params.coverage.toInteger() == 0 ? 150 : Math.round(params.coverage.toInteger() * 1.5)
    if (runtype == 'hybrid-merge-pe') {
        meta.runtype = 'hybrid'
    } else if (runtype == 'merge-pe') {
        meta.runtype = 'paired-end'
    } else if (runtype == 'merge-se') {
        meta.runtype = 'single-end'
    } else if (runtype == 'sra_accession_ont') {
        meta.runtype = 'ont'
    }
    qin = is_assembly ? 'qin=33' : 'qin=auto'
    '''
    MERGED="multiple-read-sets-merged.txt"
    mkdir -p fastqs
    mkdir -p extra

    if [ "!{runtype}" == "paired-end" ]; then
        # Paired-End Reads
        cp -L !{r1[0]} fastqs/!{meta.id}_R1.fastq.gz
        cp -L !{r2[0]} fastqs/!{meta.id}_R2.fastq.gz
        touch extra/empty.fna.gz
    elif [ "!{runtype}" == "single-end" ]; then
        # Single-End Reads
        cp -L !{r1[0]} fastqs/!{meta.id}.fastq.gz
        touch extra/empty.fna.gz
    elif [ "!{runtype}" == "ont" ]; then
        # Nanopore reads
        cp -L !{r1[0]} fastqs/!{meta.id}.fastq.gz
        touch extra/empty.fna.gz
    elif  [ "!{runtype}" == "hybrid" ]; then
        # Paired-End Reads
        cp -L !{r1[0]} fastqs/!{meta.id}_R1.fastq.gz
        cp -L !{r2[0]} fastqs/!{meta.id}_R2.fastq.gz
        cp -L !{extra} extra/!{meta.id}.fastq.gz
    elif [ "!{runtype}" == "merge-pe" ] || [ "!{runtype}" == "hybrid-merge-pe" ]; then 
        # Merge Paired-End Reads
        echo "This sample had reads merged." > ${MERGED}
        echo "R1:" >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{meta.id}_R1.fastq.gz
        echo "Merged R1:" >> ${MERGED}
        ls -l fastqs/!{meta.id}_R1.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

        echo "R2:" >> ${MERGED}
        find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
        find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{meta.id}_R2.fastq.gz
        echo "Merged R2:" >> ${MERGED}
        ls -l fastqs/!{meta.id}_R2.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

        if [ "!{runtype}" == "hybrid-merge-pe" ]; then
            cp -L !{extra} extra/!{meta.id}.fastq.gz
        else
            touch extra/empty.fna.gz
        fi
    elif [ "!{runtype}" == "merge-se" ]; then 
        # Merge Single-End Reads
        echo "This sample had reads merged." > ${MERGED}
        echo "SE:" >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{meta.id}.fastq.gz
        echo "Merged SE:" >> ${MERGED}
        ls -l fastqs/!{meta.id}.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

        touch extra/empty.fna.gz
    elif [ "!{runtype}" == "sra_accession" ] || [ "!{runtype}" == "sra_accession_ont" ]; then
        if [ "!{task.attempt}" == "!{params.max_retry}" ]; then
            echo "Unable to download !{meta.id} from both SRA and ENA !{params.max_retry} times. This may or may 
                not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                further analysis of !{meta.id} will be discontinued." | \
            sed 's/^\\s*//' > !{meta.id}-fastq-download-error.txt
            exit
        else
            # Download accession from ENA/SRA
            fastq-dl !{meta.id} !{archive} \
                --cpus !{task.cpus} \
                --outdir fastqs/ \
                --group_by_experiment \
                --is_experiment \
                --ftp_only
            touch extra/empty.fna.gz
        fi 
    elif [ "!{is_assembly}" == "true" ]; then
        if [ "!{runtype}" == "assembly_accession" ]; then
            if [ "!{task.attempt}" == "!{params.max_retry}" ]; then
                touch extra/empty.fna.gz
                echo "Unable to download !{meta.id} from NCBI Assembly !{params.max_retry} times. This may or may
                    not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                    further analysis of !{meta.id} will be discontinued." | \
                sed 's/^\\s*//' > !{meta.id}-assembly-download-error.txt
                exit
            else
                # Verify Assembly accession
                check-assembly-accession.py !{meta.id} > accession.txt 2> check-assembly-accession.txt

                if [ -s "accession.txt" ]; then
                    # Download from NCBI assembly and simulate reads
                    mkdir fasta/
                    ncbi-genome-download bacteria -o ./ -F fasta -p !{task.cpus} \
                                                -s !{section} -A accession.txt -r 50 !{no_cache}
                    find . -name "*!{meta.id}*.fna.gz" | xargs -I {} mv {} fasta/
                    rename 's/(GC[AF]_\\d+).*/$1.fna.gz/' fasta/*
                    gzip -cd fasta/!{meta.id}.fna.gz > !{meta.id}-art.fna
                    rm check-assembly-accession.txt
                else
                    mv check-assembly-accession.txt !{meta.id}-assembly-accession-error.txt
                    exit
                fi
            fi
        elif [ "!{runtype}" == "assembly" ]; then
            if [ "!{is_compressed}" == "true" ]; then
                gzip -cd !{extra} > !{meta.id}-art.fna
            else 
                cat !{extra} > !{meta.id}-art.fna
            fi
        fi

        # Simulate reads from assembly, reads are 250bp without errors
        art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov !{fcov} -ir 0 -ir2 0 -dr 0 -dr2 0 -rs !{params.sampleseed} \
                        -na -qL 33 -qU 40 -o !{meta.id}_R --id !{meta.id} -i !{meta.id}-art.fna

        mv !{meta.id}_R1.fq fastqs/!{meta.id}_R1.fastq
        mv !{meta.id}_R2.fq fastqs/!{meta.id}_R2.fastq
        pigz -p !{task.cpus} --fast fastqs/*.fastq
        cp !{meta.id}-art.fna extra/!{meta.id}.fna
        pigz -p !{task.cpus} --best extra/!{meta.id}.fna
    fi

    # Validate input FASTQs
    if [ "!{params.skip_fastq_check}" == "false" ]; then
        ERROR=0
        # Check paired-end reads have same read counts
        OPTS="--sample !{meta.id} --min_basepairs !{params.min_basepairs} --min_reads !{params.min_reads} --min_proportion !{params.min_proportion}"
        if [ -f  "fastqs/!{meta.id}_R2.fastq.gz" ]; then
            # Paired-end
            gzip -cd fastqs/!{meta.id}_R1.fastq.gz | fastq-scan > r1.json
            gzip -cd fastqs/!{meta.id}_R2.fastq.gz | fastq-scan > r2.json
            if ! reformat.sh in1=fastqs/!{meta.id}_R1.fastq.gz in2=fastqs/!{meta.id}_R2.fastq.gz !{qin} out=/dev/null 2> !{meta.id}-paired-end-error.txt; then
                ERROR=1
                echo "!{meta.id} FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." | \
                sed 's/^\\s*//' >> !{meta.id}-paired-end-error.txt
            else
                rm -f !{meta.id}-paired-end-error.txt
            fi

            if ! check-fastqs.py --fq1 r1.json --fq2 r2.json ${OPTS}; then
                ERROR=1
            fi
            rm r1.json r2.json
        else
            # Single-end
            gzip -cd fastqs/!{meta.id}.fastq.gz | fastq-scan > r1.json
            if ! check-fastqs.py --fq1 r1.json ${OPTS}; then
                ERROR=1
            fi
            rm r1.json
        fi

        # Failed validations so, let's keep them from continuing
        if [ "${ERROR}" -eq "1" ]; then
            mv fastqs/ failed-tests-fastqs/
        fi
    fi

    # Estimate Genome Size
    GENOME_SIZE_OUTPUT="!{meta.id}-genome-size.txt"
    if [ "!{meta.genome_size}" == "0" ]; then
        if [ "!{is_assembly}" == "true" ]; then
            # Use the total assembly size as the genome size
            stats.sh in=extra/!{meta.id}.fna.gz | grep All | awk '{print $5}' | sed 's/,//g' > ${GENOME_SIZE_OUTPUT}
        else
            FASTQS=""
            if [ -f  "fastqs/!{meta.id}_R2.fastq.gz" ]; then
                FASTQS="-r fastqs/!{meta.id}_R1.fastq.gz fastqs/!{meta.id}_R2.fastq.gz"
            else
                FASTQS="fastqs/!{meta.id}.fastq.gz"
            fi

            # First Pass
            mash sketch -o test -k 31 -m 3 ${FASTQS} 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${GENOME_SIZE_OUTPUT}
            rm -rf test.msh
            ESTIMATED_GENOME_SIZE=`head -n1 ${GENOME_SIZE_OUTPUT}`

            # Check if second pass is needed
            if [ ${ESTIMATED_GENOME_SIZE} -gt "!{params.max_genome_size}" ] || [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
                # Probably high coverage, try increasing number of kmer copies to 10
                M="-m 10"
                if [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
                    # Probably low coverage, try decreasing the number of kmer copies to 1
                    M="-m 1"
                fi
                mash sketch -o test -k 31 ${M} ${FASTQS} 2>&1 | \
                    grep "Estimated genome size:" | \
                    awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${GENOME_SIZE_OUTPUT}
                rm -rf test.msh
            fi
        fi

        # Check final estimate
        ESTIMATED_GENOME_SIZE=`head -n1 ${GENOME_SIZE_OUTPUT}`
        if [ ${ESTIMATED_GENOME_SIZE} -gt "!{params.max_genome_size}" ]; then
            rm ${GENOME_SIZE_OUTPUT}
            echo "!{meta.id} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) exceeds the maximum
                    allowed genome size (!{params.max_genome_size} bp). If this is unexpected, please
                    investigate !{meta.id} to determine a cause (e.g. metagenomic, contaminants, etc...).
                    Otherwise, adjust the --max_genome_size parameter to fit your need. Further analysis
                    of !{meta.id} will be discontinued." | \
            sed 's/^\\s*//' > !{meta.id}-genome-size-error.txt
        elif [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
            rm ${GENOME_SIZE_OUTPUT}
            echo "!{meta.id} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) is less than the minimum
                    allowed genome size (!{params.min_genome_size} bp). If this is unexpected, please
                    investigate !{meta.id} to determine a cause (e.g. metagenomic, contaminants, etc...).
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further analysis
                    of !{meta.id} will be discontinued." | \
            sed 's/^\\s*//' > !{meta.id}-genome-size-error.txt
        fi
    else
        # Use the genome size given by the user. (Should be >= 0)
        echo "!{meta.genome_size}" > ${GENOME_SIZE_OUTPUT}
    fi

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        art: $(echo $(art_illumina --help 2>&1) | sed 's/^.*Version //;s/ .*$//')
        fastq-dl: $(echo $(fastq-dl --version 2>&1) | sed 's/fastq-dl //')
        fastq-scan: $(echo $(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        mash: $(echo $(mash --version 2>&1))
        ncbi-genome-download: $(echo $(ncbi-genome-download --version 2>&1))
        pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    END_VERSIONS
    '''
}
