Purpose: Obtain gene counts from RNA sequencing fastqs
Programs used: HISAT2, Rsubread, samtools

Step 1: Aligning fastqs to a reference genome
  1. Download ensembl reference genome and ensemble gtf file
  2. Using HISAT2, build index files based on genome to generate SAM files and BAM files
       -Adjust python/BAM_SAM.py as necessary
       -Then, run SAM_gen.sh
  3. After above step, you should have SAM files and/or BAM files. Run csv_counts.sh to obtain gene counts using SAM files as input


Sources:

1. Kim D, Langmead B and Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nature Methods 2015
   
2. Liao Y, Smyth GK and Shi W (2019). The R package Rsubread is easier, faster,
cheaper and better for alignment and quantification of RNA sequencing reads.
Nucleic Acids Research, 47(8):e47.
http://www.ncbi.nlm.nih.gov/pubmed/30783653

3. Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics,
30(7):923-30.
http://www.ncbi.nlm.nih.gov/pubmed/24227677

4. Twelve years of SAMtools and BCFtools
Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
