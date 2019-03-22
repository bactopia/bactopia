#!/bin/bash
set -e
set -u

ariba run !{database} !{fq} !{database_name}
ariba summary !{database_name}/summary !{database_name}/report.tsv \
      --cluster_cols assembled,match,known_var,pct_id,ctg_cov,novel_var \
      --col_filter n --row_filter n
rm -rf ariba.tmp*
