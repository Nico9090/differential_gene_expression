import subprocess
import os
program_path="/path/to/programs" #Use in absence of exported path

#Output is .ht2 files which represent indices for segments of the chosen genome
#Only need to run once per genome
def build_index(organism_genome,outfile_header):
        subprocess.run([f"hisat2-build", organism_genome, outfile_header])
        return "Complete!"

#build_index("GCF_000001405.40_GRCh38.p14_genomic.fna","genome") #uncomment to run

#############################################################################

#Output is .sam which represents sequence alignment
fastq_path="/path/to/fastq/"
index_file_path="/path/to/index_files/"
def align(index_file_header,forward,reverse,output_sam): #forward and reverse are fastq files
        subprocess.run(["hisat2", "-x", index_file_header, "-1", forward, "-2", reverse, "-S", output_sam])
        return "Alignment Complete!"

align(f"{index_file_path}genome",f"{fastq_path}S12_1.fq",f"{fastq_path}S12_2.fq",
        "S12.sam")
align(f"{index_file_path}genome",f"{fastq_path}S13_1.fq",f"{fastq_path}S13_2.fq",
        "S13.sam")
align(f"{index_file_path}genome",f"{fastq_path}S14_1.fq",f"{fastq_path}S14_2.fq",
        "S14.sam")
align(f"{index_file_path}genome",f"{fastq_path}S15_1.fq",f"{fastq_path}S15_2.fq",
        "S15.sam")
align(f"{index_file_path}genome",f"{fastq_path}S16_1.fq",f"{fastq_path}S16_2.fq",
        "S16.sam")

###############################################################################################
#Need to create .bam files out of .sam to for gene counting
def create_and_sort_BAM(SAM_file,outfile,sorted_outfile):
        subprocess.run(["samtools view", "-Sb", SAM_file, ">", outfile])
        subprocess.run(["samtools sort", "-o", sorted_outfile, outfile])
        return "BAM Generated!"

SAM_list=[] #add the names of the .SAM files here once created
#for file in SAM_list:
#       file_header=file[0:file.index(".")]
#       create_and_sort_BAM(file,f"{file_header}.bam",f"{file_header}_sorted.bam")


def generate_read_counts(gtf_outfile,sorted_bam):
  subprocess.run(["stringtie", "-o ", gtf_outfile, sorted_bam])
  return None
