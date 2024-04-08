# python gen_fa.py ; ./bin/gmap cdna2.fa -g genetest2.fa -E genomic -f 3 -A

import os

from Bio import SeqIO,SeqRecord,Seq

genome = "GATGATTTTCCCGTGTGTCCCATGTAAAAAA"*1000
new_genome = SeqRecord.SeqRecord(Seq.Seq(genome),id="SOMEID",name="ultimate gene",description="hey")
cdna = SeqRecord.SeqRecord(Seq.Seq(genome[200:250:1]),id="cdna",name="kk",description="hey")
SeqIO.write(new_genome, "genetest2.fa", "fasta")
SeqIO.write(cdna, "cdna2.fa", "fasta")
os.system("./bin/gmap cdna2.fa -g genetest2.fa -E genomic -f 3 -A")