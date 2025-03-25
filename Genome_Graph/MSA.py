#Compare protein sequences with MSA
from Bio import Entrez #idea to use modules source: CHATGPT
from Bio import SeqIO
import pandas as pd
import textwrap
import subprocess
from Bio import AlignIO
import os
def search_uniprot_ids(query):
        Entrez.email = ""
        handle = Entrez.esearch(db="protein", term=query, retmax=10)  # Search for top 10 results
        record = Entrez.read(handle)
        try:
                result=record["IdList"][0]
        except IndexError:
                return None
        else:
                return result

def get_protein_sequence(accession):
    Entrez.email = ""
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return record.description, record.seq
    
if os.path.isfile("*fasta"):
    fasta_file="*.fasta"
    subprocess.run()
#for i in range(len(gene_ids)):
#        fasta_file=open((df.iloc[i]["Gene"]+".fasta"),'w')
#        try:
#                header=get_protein_sequence(gene_ids[i])[0]
#                seq=get_protein_sequence(gene_ids[i])[1]
#        except ValueError:
#                continue
#        except TypeError:
#                continue
#        else:
#                print(">",str(header),file=fasta_file)
#                print(textwrap.fill(str(seq),width=50),file=fasta_file)
#        fasta_file.close()
