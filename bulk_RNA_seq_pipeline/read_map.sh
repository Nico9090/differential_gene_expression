#!/bin/bash
hisat2 \
-x <genome_file.fa> \
{-U <*_1.fq, *_2.fq>} | \
[-S < "output.sam" >] 

