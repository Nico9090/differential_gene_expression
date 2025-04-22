#!/usr/bin/env python3
seq="GACTATATCCTAAATACCCGCACCATTACCGACACCCGTGGCCCAAGCAG"
r= len(set(seq)) #number of unique elements in seq
k=4 #length of each subsequence
max_kmers=r**k #maximum number of unique kmers you can generate from seq
start=0
all_kmers=[seq[i:i+k] for i in range(len(seq)-k+1)]
all_kmers=[kmer for kmer in all_kmers if len(kmer)==k]
print(all_kmers)
