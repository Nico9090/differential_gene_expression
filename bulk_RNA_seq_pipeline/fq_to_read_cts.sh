#!/bin/bash
#SLURM OPTIONS_____________________________________________________
#SBATCH -J "Bulk_RNA_Seq_Data_Transformation.HISAT2"
#SBATCH -o "fq_transform.out"
#SBATCH -e "fq_transform.err"
#__________________________________________________________________
#Main programs
module load python 
python align_fragments.py #script uses HISAT2 to generate SAM files
#Status check
if [ $? -ne 0 ]; then
  echo "Error from align_fragments.py execution" >&2
  exit 1
else
  echo "Process Complete!"
fi

#__________________________________________________________________
