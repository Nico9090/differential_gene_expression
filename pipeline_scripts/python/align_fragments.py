import subprocess
import os
#_________________________________________________________________________________
#Programs_________________________________________________________________________
#Step 1___________________________________________________________________________
#Example genome: GCF_000001405.40_GRCh38.p14_genomic.fna
#Example outfile_header: index, output: index.ht2, index.ht3
human_genome="GCF_000001405.40_GRCh38.p14_genomic.fna"
outfile_name="index"
def build_index(
        genome,
        outfile_header
): #run once per genome
        subprocess.run([f"hisat2-build", 
                        genome, 
                        outfile_header])
        return "Index building complete!"

build_index(human_genome,
            outfile_name)

#________________________________________________________________________________
#Step 2: Obtain SAMs by aligning FASTQs__________________________________________
fastq_path="/path/to/fastq/"
index_file_path="/path/to/index_files/" #index* files just created
def align(
        index_file_header,
        forward,
        reverse,
        output_sam
): #forward and reverse are fastq files
        subprocess.run(["hisat2", "-x", 
                        index_file_header, "-1", 
                        forward, "-2", 
                        reverse, "-S", 
                        output_sam])
        return "Alignment Complete!"

fastq_pairs=[]
fastq_pairs.append(("R1.fastq","R2.fastq"))#add pairs of fastq files
align(f"{index_file_path}{outfile_name}",
      f"{fastq_path}{fastq_pairs[0][0]}",
      f"{fastq_path}{fastq_pairs[0][1]}",
      f"{fastq_pairs[0][0][:index(".")+1]}.sam")

###############################################################################################
#Need to create .bam files out of .sam to for gene counting
def create_and_sort_BAM(SAM_file,outfile,sorted_outfile):
        with open(outfile, "wb")as out:
                 subprocess.run(["samtools", "view", "-Sb", SAM_file],stdout=out)
        with open(sorted_outfile, "wb")as sorted_out:
                subprocess.run(["samtools", "sort", "-o", sorted_outfile],stdout=sorted_out)
        return "BAM Generated!"
SAM_list=[] #add the names of the .SAM files here once created
#for file in SAM_list:
#       file_header=file[0:file.index(".")]
#       create_and_sort_BAM(file,f"{file_header}.bam",f"{file_header}_sorted.bam")


#def generate_read_counts(gtf_outfile,sorted_bam):
#  subprocess.run(["stringtie", "-o ", gtf_outfile, sorted_bam])
#  return None

def gen_counts():
        subprocess.run(["Rscript", sam_toCounts.R])
