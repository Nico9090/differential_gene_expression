import subprocess
def build_index(organism_genome,outfile_header):
  subprocess.run(["hisat2-build", organism_genome, outfile_header])
  print("Index files created!")
  subprocess.run(["ls", outfile_header + "*"])
  return None

build_index("","genome") #add path to organism genome

def align(forward,reverse,outfile_header):
  subprocess.run(["hisat2", "-x", index_file_header, "-1", fq_1, "-2", fq_2, "-S", output_sam])
  return None
align("","","output") #add forward and reverse fastq that make up the paired end reads


