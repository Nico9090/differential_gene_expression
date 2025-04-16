#!/bin/bash
#SLURM OPTIONS_____________________________________________________
#SBATCH -J "Gene Counts as CSV Generation"
#SBATCH -o "SAM_transform.out"
#SBATCH -e "SAM_transform.err"
#___________________________________________________________________
#Adjust environment if needed

#eval "$(micromamba shell hook --shell bash)"
#micromamba activate r-env
#___________________________________________________________________
#Main programs 
Rscript sam_toCounts.R
#Status check
if [ $? -ne 0 ]; then
  echo "Error from sam_toCounts.R execution" >&2
  exit 1
else
  echo "Process Complete!"
fi


#__________________________________________________________________
