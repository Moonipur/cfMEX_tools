#!/bin/bash

python ../main.py fragsize -i EE88000.hg38.frag.tsv.bgz --id EE88000
python ../main.py endmotif -i EE88000.hg38.frag.tsv.bgz --id EE88000 -t 20
