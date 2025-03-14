import subprocess
def build_index(organism_genome,outfile_header):
  subprocess.run(["hisat2-build", organism_genome, outfile_header])
  print("Index files created!")
  subprocess.run(["ls", outfile_header + "*"])
  return None

build_index("","genome") #add path to organism genome

def align(forward,reverse,outfile_header):
  subprocess.run(["hisat2", "-x", index_file_header, "-1", forward, "-2", reverse, "-S", output_sam])
  return None
align("","","output") #add forward and reverse fastq that make up the paired end reads

def create_and_sort_BAM(SAM_file,outfile,sorted_outfile):
  subprocess.run(["samtools view", "-Sb", SAM_file, ">", outfile])
  subprocess.run(["samtools sort", "-o", sorted_outfile, outfile])
  return None

create_and_sort_BAM("","","") #inpput the SAM,BAM,and,sorted BAM paths

def generate_read_counts(gtf_outfile,sorted_bam):
  subprocess.run(["stringtie", "-o ", gtf_outfile, sorted_bam])
  return None
