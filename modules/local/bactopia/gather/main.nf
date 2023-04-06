// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES      = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options        = initOptions(params.options ? params.options : [:], 'gather')
options.ignore = [".fastq.gz", ".fna.gz"]
options.btype  = options.btype ?: "main"
conda_tools    = "bioconda::bactopia-gather=1.0.1"
conda_name     = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env      = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process GATHER {
    tag "${meta.id}"
    label "process_low"
    maxForks params.max_downloads
    maxRetries params.max_retry

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-gather:1.0.1--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-gather:1.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(r1, stageAs: '*???-r1'), path(r2, stageAs: '*???-r2'), path(extra)

    output:
    tuple val(meta), path("fastqs/${prefix}*.fastq.gz"), path("extra/*.gz"), emit: raw_fastq, optional: true
    tuple val(meta), path("fastqs/${prefix}*.fastq.gz"), emit: fastq_only, optional: true
    tuple val(meta), path("${prefix}-meta.tsv")        , emit: tsv
    path "*.{log,err}"         , emit: logs, optional: true
    path ".command.*"          , emit: nf_logs
    path "versions.yml"        , emit: versions
    path "*-{error,merged}.txt", optional: true

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    meta.original_runtype = meta.runtype
    runtype = meta.original_runtype
    is_assembly = runtype.startsWith('assembly') ? true : false
    is_compressed = extra ? (extra.getName().endsWith('gz') ? true : false) : false
    no_cache = params.no_cache ? '-N' : ''
    archive = params.use_ena ? (task.attempt >= 4 ? "SRA" : "ENA") : "SRA"
    section = runtype == 'assembly_accession' ? (prefix.startsWith('GCF') ? 'refseq' : 'genbank') : null
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
    """
    MERGED="multiple-read-sets-merged.txt"
    mkdir -p fastqs
    mkdir -p extra

    if [ "${runtype}" == "paired-end" ]; then
        # Paired-End Reads
        cp -L ${r1[0]} fastqs/${prefix}_R1.fastq.gz
        cp -L ${r2[0]} fastqs/${prefix}_R2.fastq.gz
        touch extra/empty.fna.gz
    elif [ "${runtype}" == "single-end" ]; then
        # Single-End Reads
        cp -L ${r1[0]} fastqs/${prefix}.fastq.gz
        touch extra/empty.fna.gz
    elif [ "${runtype}" == "ont" ]; then
        # Nanopore reads
        cp -L ${r1[0]} fastqs/${prefix}.fastq.gz
        touch extra/empty.fna.gz
    elif  [ "${runtype}" == "hybrid" ] || [ "${runtype}" == "short_polish" ]; then 
        # Paired-End Reads
        cp -L ${r1[0]} fastqs/${prefix}_R1.fastq.gz
        cp -L ${r2[0]} fastqs/${prefix}_R2.fastq.gz
        cp -L ${extra} extra/${prefix}.fastq.gz
    elif [ "${runtype}" == "merge-pe" ] || [ "${runtype}" == "hybrid-merge-pe" ]; then 
        # Merge Paired-End Reads
        echo "This sample had reads merged." > \${MERGED}
        echo "R1:" >> \${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print \$5"\t"\$9}' >> \${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/${prefix}_R1.fastq.gz
        echo "Merged R1:" >> \${MERGED}
        ls -l fastqs/${prefix}_R1.fastq.gz | awk '{print \$5"\t"\$9}' >> \${MERGED}

        echo "R2:" >> \${MERGED}
        find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print \$5"\t"\$9}' >> \${MERGED}
        find -name "*r2" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/${prefix}_R2.fastq.gz
        echo "Merged R2:" >> \${MERGED}
        ls -l fastqs/${prefix}_R2.fastq.gz | awk '{print \$5"\t"\$9}' >> \${MERGED}

        if [ "${runtype}" == "hybrid-merge-pe" ]; then
            cp -L ${extra} extra/${prefix}.fastq.gz
        else
            touch extra/empty.fna.gz
        fi
    elif [ "${runtype}" == "merge-se" ]; then 
        # Merge Single-End Reads
        echo "This sample had reads merged." > \${MERGED}
        echo "SE:" >> \${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} ls -l {} | awk '{print \$5"\t"\$9}' >> \${MERGED}
        find -name "*r1" | sort | xargs -I {} readlink {} | xargs -I {} cat {} > fastqs/${prefix}.fastq.gz
        echo "Merged SE:" >> \${MERGED}
        ls -l fastqs/${prefix}.fastq.gz | awk '{print \$5"\t"\$9}' >> \${MERGED}

        touch extra/empty.fna.gz
    elif [ "${runtype}" == "sra_accession" ] || [ "${runtype}" == "sra_accession_ont" ]; then
        if [ "${task.attempt}" == "${params.max_retry}" ]; then
            echo "Unable to download ${prefix} from both SRA and ENA ${params.max_retry} times. This may or may 
                not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                further analysis of ${prefix} will be discontinued." | \\
            sed 's/^\\s*//' > ${prefix}-fastq-download-error.txt
            exit
        else
            # Download accession from ENA/SRA
            fastq-dl \\
                --accession ${prefix} \\
                --provider ${archive} \\
                --cpus ${task.cpus} \\
                --outdir fastqs/ \\
                --prefix ${prefix} \\
                --group-by-experiment
            touch extra/empty.fna.gz
        fi 
    elif [ "${is_assembly}" == "true" ]; then
        if [ "${runtype}" == "assembly_accession" ]; then
            if [ "${task.attempt}" == "${params.max_retry}" ]; then
                touch extra/empty.fna.gz
                echo "Unable to download ${prefix} from NCBI Assembly ${params.max_retry} times. This may or may
                    not be a temporary connection issue. Rather than stop the whole Bactopia run, 
                    further analysis of ${prefix} will be discontinued." | \\
                sed 's/^\\s*//' > ${prefix}-assembly-download-error.txt
                exit
            else
                # Verify Assembly accession
                check-assembly-accession.py ${prefix} > accession.txt 2> check-assembly-accession.txt

                if [ -s "accession.txt" ]; then
                    # Download from NCBI assembly and simulate reads
                    mkdir fasta/
                    ncbi-genome-download bacteria -o ./ -F fasta -p ${task.cpus} \\
                                                -s ${section} -A accession.txt -r 50 ${no_cache}
                    find . -name "*${prefix}*.fna.gz" | xargs -I {} mv {} fasta/
                    rename 's/(GC[AF]_\\d+).*/\$1.fna.gz/' fasta/*
                    gzip -cd fasta/${prefix}.fna.gz > ${prefix}-art.fna
                    rm check-assembly-accession.txt
                else
                    mv check-assembly-accession.txt ${prefix}-assembly-accession-error.txt
                    exit
                fi
            fi
        elif [ "${runtype}" == "assembly" ]; then
            if [ "${is_compressed}" == "true" ]; then
                gzip -cd ${extra} > ${prefix}-art.fna
            else 
                cat ${extra} > ${prefix}-art.fna
            fi
        fi

        # Simulate reads from assembly, reads are 250bp without errors
        art_illumina -p -ss MSv3 -l 250 -m 400 -s 30 --fcov ${fcov} -ir 0 -ir2 0 -dr 0 -dr2 0 -rs ${params.sampleseed}\
                        -na -qL 33 -qU 40 -o ${prefix}_R --id ${prefix} -i ${prefix}-art.fna

        mv ${prefix}_R1.fq fastqs/${prefix}_R1.fastq
        mv ${prefix}_R2.fq fastqs/${prefix}_R2.fastq
        pigz -p ${task.cpus} --fast fastqs/*.fastq
        cp ${prefix}-art.fna extra/${prefix}.fna
        pigz -p ${task.cpus} --best extra/${prefix}.fna
    fi

    # Validate input FASTQs
    if [ "${params.skip_fastq_check}" == "false" ]; then
        ERROR=0
        # Check paired-end reads have same read counts
        OPTS="--sample ${prefix} --min_basepairs ${params.min_basepairs} --min_reads ${params.min_reads} --min_proportion ${params.min_proportion} --runtype ${runtype}"
        if [ -f  "fastqs/${prefix}_R2.fastq.gz" ]; then
            # Paired-end
            gzip -cd fastqs/${prefix}_R1.fastq.gz | fastq-scan > r1.json
            gzip -cd fastqs/${prefix}_R2.fastq.gz | fastq-scan > r2.json
            if ! reformat.sh in1=fastqs/${prefix}_R1.fastq.gz in2=fastqs/${prefix}_R2.fastq.gz ${qin} out=/dev/null 2> ${prefix}-paired-end-error.txt; then
                ERROR=1
                echo "${prefix} FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." | \\
                sed 's/^\\s*//' >> ${prefix}-paired-end-error.txt
            else
                rm -f ${prefix}-paired-end-error.txt
            fi

            if ! check-fastqs.py --fq1 r1.json --fq2 r2.json \${OPTS}; then
                ERROR=1
            fi
            rm r1.json r2.json
        else
            # Single-end
            gzip -cd fastqs/${prefix}.fastq.gz | fastq-scan > r1.json
            if ! check-fastqs.py --fq1 r1.json \${OPTS}; then
                ERROR=1
            fi
            rm r1.json
        fi

        # Failed validations so, let's keep them from continuing
        if [ "\${ERROR}" -eq "1" ]; then
            mv fastqs/ failed-tests-fastqs/
        fi
    fi

    # Dump meta values to a TSV
    echo "sample<TAB>runtype<TAB>original_runtype<TAB>species<TAB>genome_size" | sed 's/<TAB>/\t/g' > ${prefix}-meta.tsv
    echo "${meta.id}<TAB>${meta.runtype}<TAB>${meta.original_runtype}<TAB>${meta.species}<TAB>${meta.genome_size}" | sed 's/<TAB>/\t/g' >> ${prefix}-meta.tsv

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        art: \$(echo \$(art_illumina --help 2>&1) | sed 's/^.*Version //;s/ .*\$//')
        fastq-dl: \$(echo \$(fastq-dl --version 2>&1) | sed 's/fastq-dl //')
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        ncbi-genome-download: \$(echo \$(ncbi-genome-download --version 2>&1))
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/pigz //')
    END_VERSIONS
    """
}
