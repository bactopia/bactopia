#!/bin/bash
set -e
set -u

# Get Mash distance
mash dist -t !{refseq_sketch} !{sample_sketch} | sort -k 3,3 > mash-dist.txt

# Pick genomes to download
!{baseDir}/bin/select-references.py mash-dist.txt !{max_references} !{tie_break} > accession-list.txt

# Download genomes
ncbi-genome-download bacteria -l complete -o ./ -F genbank -p !{task.cpus} -A accession-list.txt

# Split Genbank files containing plasmids
mkdir -p refseq/split
ls refseq/bacteria/ | \
    xargs -I {} sh -c 'zcat refseq/bacteria/{}/{}*.gbff.gz | csplit -z --quiet --prefix=refseq/split/{}- - "/^\/\//+1" "{*}"'

# Copy the largest segment for each split Genbank (assumes completed genomes)
mkdir genbank
ls refseq/bacteria/ | \
    xargs -I {} sh -c 'ls -l -S refseq/split/{}-* | head -n 1 | awk "{print \$9\" genbank/{}.gbk\"}"' | \
    xargs -I {} sh -c 'cp {}'
