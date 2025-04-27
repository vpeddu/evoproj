#!/bin/bash

    ls -lah
    base=$(basename "LTR78_human.bed")
    IFS='_' read -r -a parts <<< "$base"
    repeat="${parts[0]}"
    species="${parts[-1]%.bed}"

python3 /media/vikas/fast/data/TEBAG/publication/scripts/evoproj/bin/permutation_test.py     "LTR78_human.bed"     .     chm13v2.0.fa.fai     "${repeat}"     "${species}"     "H3K4me3"     10000
