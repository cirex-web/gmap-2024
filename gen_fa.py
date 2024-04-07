# python gen_fa.py ; ./bin/gmap cdna2.fa -g genetest2.fa -E genomic -f 3 -A


from Bio import SeqIO,SeqRecord,Seq

new_genome = SeqRecord.SeqRecord(Seq.Seq("AAAAAAAAAA"*100),id="SOMEID",name="ultimate gene",description="hey");
cdna = SeqRecord.SeqRecord(Seq.Seq("GGGGGG"+"AAAAAAAAAAAA"*5+"GGGG"),id="cdna",name="kk",description="hey");
SeqIO.write(new_genome, "genetest2.fa", "fasta")
SeqIO.write(cdna, "cdna2.fa", "fasta")
