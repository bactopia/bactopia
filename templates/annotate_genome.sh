#!/bin/bash
set -e
set -u

gunzip -f !{fasta}
prokka --cpus !{cpus} --outdir annotation --force !{proteins} \
    --prefix !{sample} --locustag !{sample} --centre !{params.centre} \
    --addgenes --mincontiglen 500 !{gunzip_fasta}

find annotation/ -type f -not -name "*.txt" -and -not -name "*.log*" | \
    xargs -I {} pigz -n -p !{cpus} {}
