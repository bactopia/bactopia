#!/bin/bash
set -e
set -u

if [ "!{params.dry_run}" == "true" ]; then
    mkdir -p insertion-sequences
    touch insertion-sequences/!{sample}-!{insertion_name}.log
else
      if [[ !{params.compress} == "true" ]]; then
          gunzip -f !{genbank}
      fi

      ismap --reads !{sample}_R*.fastq.gz \
            --queries !{insertion_fasta} \
            --reference !{gunzip_genbank} \
            --log !{sample}-!{insertion_name} \
            --min_clip !{params.min_clip} \
            --max_clip !{params.max_clip} \
            --cutoff !{params.cutoff} \
            --novel_gap_size !{params.novel_gap_size} \
            --min_range !{params.min_range} \
            --max_range !{params.max_range} \
            --merging !{params.merging} \
            --T !{params.ismap_minqual} \
            --t !{task.cpus} !{all}

      mv !{sample} insertion-sequences
      mv !{sample}-!{insertion_name}.log insertion-sequences
fi
