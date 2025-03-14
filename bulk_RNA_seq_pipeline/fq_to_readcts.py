import subprocess
genome_file="" #path to genome file
index_file_header="genome" #output of building index files with hisat2
subprocess.run(["hisat2-build", genome_file, index_file_header]) #create index for read alignment
print("Index files created!")
subprocess.run(["ls",index_file_header*])
