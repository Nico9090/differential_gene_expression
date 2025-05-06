#!/bin/bash
module load fasterq-dump
INDIR=""
OUTDIR=""
cd $INDIR
SRR_file=""
fasterq-dump --outdir $OUTDIR --split-3 $SRR_file
