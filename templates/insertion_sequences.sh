#!/bin/bash
set -e
set -u

gunzip -f !{genbank}
ismap --reads !{sample}_R*.fastq.gz --queries !{insertion_fasta} \
      --reference !{gunzip_genbank} --log !{sample}-!{insertion_name} \
      --t !{task.cpus}
mv !{sample} insertion-sequences
mv !{sample}-!{insertion_name}.log insertion-sequences
