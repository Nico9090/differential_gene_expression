#!/bin/bash
#SLURM OPTIONS_____________________________________________________
#SBATCH -J "SAM_File_Creation_Using_HISAT2"
#SBATCH -o "fastq_transform.out"
#SBATCH -e "fastq_transform.err"
#___________________________________________________________________
#Adjust environment if needed

#eval "$(micromamba shell hook --shell bash)"
#micromamba activate bulk_RNA_seq 

#___________________________________________________________________
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
