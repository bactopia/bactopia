nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "gather_samples"

process GATHER_SAMPLES {
    /* Gather up input FASTQs for analysis. */
    publishDir "${params.outdir}/${sample}",
        mode: params.publish_mode,
        overwrite: params.overwrite,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME, ignore: [".fastq.gz", ".fna.gz"]) }

    tag "${sample}"
    label "max_cpus"
    label PROCESS_NAME

    input:
    tuple val(sample), val(sample_type), val(genome_size), file(r1: '*???-r1'), file(r2: '*???-r2'), path(extra)

    output:
    tuple val(sample), val(final_sample_type), path("fastqs/${sample}*.fastq.gz"),
          path("extra/*.gz"), path("${sample}-genome-size.txt"), emit: raw_fastq, optional: true
    path "*.std{out,err}.txt", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version
    path "*-{error,merged}.txt", optional:true

    shell:
    bactopia_version = workflow.manifest.version
    nextflow_version = nextflow.version
    is_assembly = sample_type.startsWith('assembly') ? true : false
    is_compressed = extra ? (extra.getName().endsWith('gz') ? true : false) : false
    no_cache = params.no_cache ? '-N' : ''
    archive = params.use_ena ? (task.attempt >= 4 ? "SRA" : "ENA") : "SRA"
    section = sample_type == 'assembly-accession' ? (sample.startsWith('GCF') ? 'refseq' : 'genbank') : null
    fcov = params.coverage.toInteger() == 0 ? 150 : Math.round(params.coverage.toInteger() * 1.5)
    final_sample_type = sample_type
    if (sample_type == 'hybrid-merge-pe') {
        final_sample_type = 'hybrid'
    } else if (sample_type == 'merge-pe') {
        final_sample_type = 'paired-end'
    } else if (sample_type == 'merge-se') {
        final_sample_type = 'single-end'
    }
    qin = is_assembly ? 'qin=33' : 'qin=auto'
    '''
    MERGED="multiple-read-sets-merged.txt"
    mkdir -p fastqs
    mkdir -p extra

    if [ "!{sample_type}" == "paired-end" ]; then
        # Paired-End Reads
        ln -s `readlink !{r1[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{r2[0]}` fastqs/!{sample}_R2.fastq.gz
        touch extra/empty.fna.gz
    elif [ "!{sample_type}" == "single-end" ]; then
        # Single-End Reads
        ln -s `readlink !{r1[0]}` fastqs/!{sample}.fastq.gz
        touch extra/empty.fna.gz
    elif  [ "!{sample_type}" == "hybrid" ]; then
        # Paired-End Reads
        ln -s `readlink !{r1[0]}` fastqs/!{sample}_R1.fastq.gz
        ln -s `readlink !{r2[0]}` fastqs/!{sample}_R2.fastq.gz
        ln -s `readlink !{extra}` extra/!{sample}.fastq.gz
    elif [ "!{sample_type}" == "merge-pe" ] || [ "!{sample_type}" == "hybrid-merge-pe" ]; then 
        # Merge Paired-End Reads
        echo "This sample had reads merged." > ${MERGED}
        echo "R1:" >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}_R1.fastq.gz
        echo "Merged R1:" >> ${MERGED}
        ls -l fastqs/!{sample}_R1.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

        echo "R2:" >> ${MERGED}
        find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
        find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}_R2.fastq.gz
        echo "Merged R2:" >> ${MERGED}
        ls -l fastqs/!{sample}_R2.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

        if [ "!{sample_type}" == "hybrid-merge-pe" ]; then
            ln -s `readlink !{extra}` extra/!{sample}.fastq.gz
        else
            touch extra/empty.fna.gz
        fi
    elif [ "!{sample_type}" == "merge-se" ]; then 
        # Merge Single-End Reads
        echo "This sample had reads merged." > ${MERGED}
        echo "SE:" >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print $5"\t"$9}' >> ${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/!{sample}.fastq.gz
        echo "Merged SE:" >> ${MERGED}
        ls -l fastqs/!{sample}.fastq.gz | awk '{print $5"\t"$9}' >> ${MERGED}

        touch extra/empty.fna.gz
    elif [ "!{sample_type}" == "sra-accession" ]; then
        # fastq-dl Version
        fastq-dl --version > fastq-dl.version.txt 2>&1

        if [ "!{task.attempt}" == "!{params.max_retry}" ]; then
            echo "Unable to download !{sample} from both SRA and ENA !{params.max_retry} times. This may or may 
                not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                further analysis of !{sample} will be discontinued." | \
            sed 's/^\\s*//' > !{sample}-fastq-download-error.txt
            exit
        else
            # Download accession from ENA/SRA
            fastq-dl !{sample} !{archive} \
                --cpus !{task.cpus} \
                --outdir fastqs/ \
                --group_by_experiment \
                --is_experiment \
                --ftp_only  > fastq-dl.stdout.txt 2> fastq-dl.stderr.txt
            touch extra/empty.fna.gz
        fi 
    elif [ "!{is_assembly}" == "true" ]; then
        if [ "!{sample_type}" == "assembly-accession" ]; then
            # ncbi-genome-download Version
            ncbi-genome-download --version > ncbi-genome-download.version.txt 2>&1

            if [ "!{task.attempt}" == "!{params.max_retry}" ]; then
                touch extra/empty.fna.gz
                echo "Unable to download !{sample} from NCBI Assembly !{params.max_retry} times. This may or may
                    not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                    further analysis of !{sample} will be discontinued." | \
                sed 's/^\\s*//' > !{sample}-assembly-download-error.txt
                exit
            else
                # Verify Assembly accession
                check-assembly-accession.py !{sample} > accession.txt 2> check-assembly-accession.stderr.txt

                if [ -s "accession.txt" ]; then
                    # Download from NCBI assembly and simulate reads
                    mkdir fasta/
                    ncbi-genome-download bacteria -o ./ -F fasta -p !{task.cpus} \
                                                -s !{section} -A accession.txt -r 50 !{no_cache} > ncbi-genome-download.stdout.txt 2> ncbi-genome-download.stderr.txt
                    find . -name "*!{sample}*.fna.gz" | xargs -I {} mv {} fasta/
                    rename 's/(GC[AF]_\\d+).*/$1.fna.gz/' fasta/*
                    gzip -cd fasta/!{sample}.fna.gz > !{sample}-art.fna
                else
                    cp check-assembly-accession.stderr.txt !{sample}-assembly-accession-error.txt
                    exit
                fi
            fi
        elif [ "!{sample_type}" == "assembly" ]; then
            if [ "!{is_compressed}" == "true" ]; then
                gzip -cd !{extra} > !{sample}-art.fna
            else 
                cat !{extra} > !{sample}-art.fna
            fi
        fi

        # Simulate reads from assembly, reads are 250bp without errors
        art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov !{fcov} -ir 0 -ir2 0 -dr 0 -dr2 0 -rs !{params.sampleseed} \
                        -na -qL 33 -qU 40 -o !{sample}_R --id !{sample} -i !{sample}-art.fna > art.stdout.txt 2> art.stderr.txt

        mv !{sample}_R1.fq fastqs/!{sample}_R1.fastq
        mv !{sample}_R2.fq fastqs/!{sample}_R2.fastq
        pigz -p !{task.cpus} --fast fastqs/*.fastq
        cp !{sample}-art.fna extra/!{sample}.fna
        pigz -p !{task.cpus} --best extra/!{sample}.fna

        # ART Version
        art_illumina --help | head -n 6 | tail -n 5 > art.version.txt 2>&1
    fi

    # Validate input FASTQs
    if [ "!{params.skip_fastq_check}" == "false" ]; then
        ERROR=0
        # Not completely sure about the inputs, so make sure they meet minimum requirements
        fastq-scan -v > fastq-scan.version.txt 2>&1

        # Check paired-end reads have same read counts
        OPTS="--sample !{sample} --min_basepairs !{params.min_basepairs} --min_reads !{params.min_reads} --min_proportion !{params.min_proportion}"
        if [ -f  "fastqs/!{sample}_R2.fastq.gz" ]; then
            # Paired-end
            gzip -cd fastqs/!{sample}_R1.fastq.gz | fastq-scan > r1.json
            gzip -cd fastqs/!{sample}_R2.fastq.gz | fastq-scan > r2.json
            if ! reformat.sh in1=fastqs/!{sample}_R1.fastq.gz in2=fastqs/!{sample}_R2.fastq.gz !{qin} out=/dev/null 2> !{sample}-paired-end-error.txt; then
                ERROR=1
                echo "!{sample} FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." | \
                sed 's/^\\s*//' >> !{sample}-paired-end-error.txt
            else
                rm -f !{sample}-paired-end-error.txt
            fi

            if ! check-fastqs.py --fq1 r1.json --fq2 r2.json ${OPTS}; then
                ERROR=1
            fi
            rm r1.json r2.json
        else
            # Single-end
            gzip -cd fastqs/!{sample}.fastq.gz | fastq-scan > r1.json
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
    GENOME_SIZE_OUTPUT="!{sample}-genome-size.txt"
    if [ "!{genome_size}" == "0" ]; then
        if [ "!{is_assembly}" == "true" ]; then
            # Use the total assembly size as the genome size
            stats.sh in=extra/!{sample}.fna.gz | grep All | awk '{print $5}' | sed 's/,//g' > ${GENOME_SIZE_OUTPUT}
        else
            FASTQS=""
            if [ -f  "fastqs/!{sample}_R2.fastq.gz" ]; then
                FASTQS="-r fastqs/!{sample}_R1.fastq.gz fastqs/!{sample}_R2.fastq.gz"
            else
                FASTQS="fastqs/!{sample}.fastq.gz"
            fi
            # Use mash to estimate the genome size, if a genome size cannot be estimated set the genome size to 0
            mash --version > mash.version.txt 2>&1

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
            echo "!{sample} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) exceeds the maximum
                    allowed genome size (!{params.max_genome_size} bp). If this is unexpected, please
                    investigate !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...).
                    Otherwise, adjust the --max_genome_size parameter to fit your need. Further analysis
                    of !{sample} will be discontinued." | \
            sed 's/^\\s*//' > !{sample}-genome-size-error.txt
        elif [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
            rm ${GENOME_SIZE_OUTPUT}
            echo "!{sample} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) is less than the minimum
                    allowed genome size (!{params.min_genome_size} bp). If this is unexpected, please
                    investigate !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...).
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further analysis
                    of !{sample} will be discontinued." | \
            sed 's/^\\s*//' > !{sample}-genome-size-error.txt
        fi
    else
        # Use the genome size given by the user. (Should be >= 0)
        echo "!{genome_size}" > ${GENOME_SIZE_OUTPUT}
    fi

    # Capture versions
    echo "bactopia !{bactopia_version}" > bactopia.version.txt
    echo "nextflow !{nextflow_version}" > nextflow.version.txt
    '''

    stub:
    final_sample_type = 'single-end'
    """
    mkdir fastqs
    mkdir extra
    touch ${sample}-error.txt
    touch fastqs/${sample}.fastq.gz
    touch extra/${sample}.gz
    touch bactopia.versions
    touch multiple-read-sets-merged.txt
    """
}
