#!/bin/bash
set -e
set -u
LOG_DIR="!{task.process}"
mkdir -p ${LOG_DIR}
echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions

# Verify AWS files were staged
if [[ ! -L "!{fq[0]}" ]]; then
    if [ "!{single_end}" == "true" ]; then
        check_staging.py --fq1 !{fq[0]} --extra !{sample_sketch} --is_single
    else
        check_staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{sample_sketch}
    fi
fi

# Get Mash distance
echo "# Mash Version" >> ${LOG_DIR}/!{task.process}.versions
mash --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
mash dist -t !{sample_sketch} !{refseq_sketch} | grep -v "query" | sort -k 2,2 > distances.txt

# Pick genomes to download
printf "accession\tdistance\tlatest_accession\tupdated\n" > mash-dist.txt
select-references.py distances.txt !{total} !{tie_break} >> mash-dist.txt

# Pick only latest accessions
grep -v distance mash-dist.txt | cut -f3 > download-list.txt

# Download genomes
echo "# ncbi-genome-download Version" >> ${LOG_DIR}/!{task.process}.versions
ncbi-genome-download --version >> ${LOG_DIR}/!{task.process}.versions 2>&1
ncbi-genome-download bacteria -l complete -o ./ -F genbank -p !{task.cpus} -A download-list.txt -r 50 !{no_cache} > ${LOG_DIR}/ncbi-genome-download.out 2> ${LOG_DIR}/ncbi-genome-download.err

# Move and uncompress genomes
mkdir genbank_temp
find refseq -name "*.gbff.gz" | xargs -I {} mv {} genbank_temp/
rename 's/(GC[AF]_\d+).*/$1/' genbank_temp/*
mkdir genbank
ls genbank_temp/ | xargs -I {} sh -c 'gzip -cd genbank_temp/{} > genbank/!{sample}-{}.gbk'
rm -rf genbank_temp

if [ "!{params.keep_all_files}" == "false" ]; then
    # Remove intermediate GenBank files
    rm -rf refseq/
fi

if [ "!{params.skip_logs}" == "false" ]; then 
    cp .command.err ${LOG_DIR}/!{task.process}.err
    cp .command.out ${LOG_DIR}/!{task.process}.out
    cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
    cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
else
    rm -rf ${LOG_DIR}/
fi
