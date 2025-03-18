#!/bin/bash
#SBATCH -J "Bulk_RNA_Seq_Data_Transformation.HISAT2_STRINGTIE"
#SBATCH -o "fq_transform.out"
#SBATCH -e "fq_transform.err"
#SBATCH --partition=long
#SBATCH --array=0-17
#SBATCH --nodes=1
#SBATCH --cpus-per-task=70
#SBATCH --mem=300G

module load python
python align_fragments.py
#/path/to/programs/hisat2/hisat2-build /path/to/programs/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna genome
echo "Process Complete!"
