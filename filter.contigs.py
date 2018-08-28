import csv
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('infasta', help='sorted bam file with mapped reads')
ap.add_argument('ids',help='list of ids which were mapped to the human reference')
ap.add_argument('out', help='file to save the number of reads per genome category')
args = ap.parse_args()

ids=set()

file=open(args.ids)
reader=csv.reader(file)
for line in reader:
	ids.add(line[0])

print "Microbial contigs mapped to human", ids


from Bio import SeqIO



fasta_sequences = SeqIO.parse(open(args.infasta),'fasta')
out_file=open(args.out,"w")
k=0

for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name not  in ids:
            k+=1
            out_file.write(">"+name)
            out_file.write("\n")
            out_file.write(sequence)
            out_file.write("\n")






out_file.close()

print k, "microbial contigs are non-human"
print "Success!",


