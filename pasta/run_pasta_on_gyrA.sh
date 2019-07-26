#!/bin/bash

# Set path to SEPP program
PASTA="/opt/pasta-code/pasta/run_pasta.py"

# Set path to RDP 2016 Reference Package
RAW_SEQUENCES="./gyrA_raw.fasta"

# Decide on 
JOB_NAME="gryA_pasta"

# Run PASTA
python $PASTA -i $RAW_SEQUENCES \
             --temporaries "./tmp/" \
             -d Protein \
             -j $JOB_NAME
 