wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
gunzip *.gz #unzip fastq files as well if they are zipped

#make genome index file
STAR --runThreadN 64 --runMode genomeGenerate --genomeDir /data/genomes/h38/STAR/ --genomeFastaFiles /data/genomes/h38/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /data/genomes/h38/STAR/Homo_sapiens.GRCh38.99.gtf

#source: https://github.com/erilu/bulk-rnaseq-analysis/blob/master/README.md
