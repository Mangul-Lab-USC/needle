import argparse
import re
import sys
import csv


ap = argparse.ArgumentParser()
ap.add_argument('input_blast', help='--')
ap.add_argument('input_BWA', help='--')
ap.add_argument('out', help='---')
args = ap.parse_args()




file = open(args.input_blast,"r")



#contig,id,name,identity,alignment_length,contig_length,adjusted_identity
#k21_2,AB261990.1,Lymphocytic_choriomeningitis_virus_genomic_RNA_segment_S_nearly_,100.0,916.0,916.0,100.0



print ("Open ",args.input_blast)

file=open(args.input_blast,"r")
reader=csv.reader(file)
next(reader,None)

dict_blast={}
kmers=set()

for line in reader:
    kmer_name=line[0]
    ref_name_short=line[2].split("_")[0]
    identity=float(line[6])
    if kmer_name not in kmers: #kmer occurs first time => the best hit
        kmers.add(kmer_name)
        dict_blast[kmer_name] = line


file.close()






print ("Open ",args.input_BWA)
#contig,id,name,number_mismatches,alignment_length,contig_length,adjusted_identity
#k21_2,gb|AB261990|Strain|M2|Description|Lymphocytic,gb|AB261990|Strain|M2|Description|Lymphocytic,0,916,916,100.0

file=open(args.input_BWA,"r")
reader=csv.reader(file)
next(reader,None)

dict_BWA={}
kmers=set()

for line in reader:
    kmer_name=line[0]
    ref_name_short=line[1].split("|")[1]
    identity=float(line[6])
    if kmer_name not in kmers: #kmer occurs first time => the best hit
        kmers.add(kmer_name)
        dict_BWA[kmer_name] = line


file.close()


print (dict_blast)
print (dict_BWA)


print set(dict_BWA)-set(dict_blast)



fileOut=open(args.out,"w")
fileOut.write("contig,id,name,number_mismatches,alignment_length,contig_length,adjusted_identity\n")


for i in set(dict_BWA)-set(dict_blast):
    str=','.join(dict_BWA[i])
    fileOut.write(str)
    fileOut.write("\n")








