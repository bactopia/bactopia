#!/bin/bash
set -e
set -u

if [ "!{params.dry_run}" == "true" ]; then
    touch mash-dist.txt
else
    # Get Mash distance
    mash dist -t !{sample_sketch} !{refseq_sketch} | sort -k 2,2 > distances.txt

    # Pick genomes to download
    printf "accession\tdistance\tlatest_accession\tupdated\n" > mash-dist.txt
    select-references.py distances.txt !{total} !{tie_break} >> mash-dist.txt

    # Pick only latest accessions
    grep -v distance mash-dist.txt | cut -f3 > download-list.txt

    # Download genomes
    ncbi-genome-download bacteria -l complete -o ./ -F genbank -p !{task.cpus} -A download-list.txt -r 50

    # Split Genbank files containing plasmids
    mkdir -p refseq/split
    ls refseq/bacteria/ | \
        xargs -I {} sh -c 'zcat refseq/bacteria/{}/{}*.gbff.gz | csplit -z --quiet --prefix=refseq/split/{}- - "/^\/\//+1" "{*}"'

    # Copy the largest segment for each split Genbank (assumes completed genomes)
    mkdir genbank
    ls refseq/bacteria/ | \
        xargs -I {} sh -c 'ls -l -S refseq/split/{}-* | head -n 1 | awk "{print \$9\" genbank/!{sample}-{}.gbk\"}"' | \
        xargs -I {} sh -c 'cp {}'

    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove intermediate GenBank files
        rm -rf refseq/
    fi
fi
