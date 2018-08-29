import argparse
import re
import sys
import csv
import gzip
import os


#comapre indetify of 10 hits from BLAST
#['k21_2', 'AB261990.1', 'Lymphocytic_choriomeningitis_virus_genomic_RNA_segment_S_nearly_', '100.0', '916.0', '916.0', '100.0']
def compare(l):



    best_genus=l[0][2].split("_")[0]
    best_identity=identity=float(l[0][6])

    for i in l[1:]:
        identity=float(i[6])
        genus_blast = i[2].split("_")[0]
        id_blast =i[1].split(".")[0]
        if genus_blast!=best_genus:
            if abs(identity-best_identity)<2.0:
		print "compare--",l,identity,best_identity,identity-best_identity
                return False # we have found hit from other genus with within 3% similarity

    return True # we have NOT found hit from other genus with within 3% similarity





ap = argparse.ArgumentParser()
ap.add_argument('input_blast', help='--')
ap.add_argument('input_BWA', help='--')
ap.add_argument('out', help='---')
#ap.add_argument("-o",help="[virus, fungi or protozoa]. By default virus", default='virus',type=str)
args = ap.parse_args()





file = open(args.input_blast,"r")



#contig,id,name,identity,alignment_length,contig_length,adjusted_identity
#k21_2,AB261990.1,Lymphocytic_choriomeningitis_virus_genomic_RNA_segment_S_nearly_,100.0,916.0,916.0,100.0



print ("Open ",args.input_blast)

file=open(args.input_blast,"r")
reader=csv.reader(file)
next(reader,None)

dict_blast={}
dict_blast2={} # to store secondary alligments
kmers=set()

for line in reader:
    kmer_name=line[0]
    ref_name_short=line[2].split("_")[0]
    identity=float(line[6])
    if kmer_name not in kmers: #kmer occurs first time => the best hit
        kmers.add(kmer_name)
        dict_blast[kmer_name] = line
        dict_blast2[kmer_name] = []
        dict_blast2[kmer_name].append(line)
    else:
        dict_blast2[kmer_name].append(line)


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
    ref_id=line[1]
    identity=float(line[6])
    if kmer_name not in kmers: #kmer occurs first time => the best hit
        kmers.add(kmer_name)
        dict_BWA[kmer_name] = line



file.close()





print "Contigs which are only mapped to custom microbial database and is not mapped to TREE OF LIFE FROM BLAST"
print set(dict_BWA)-set(dict_blast)



fileOut_blast=open(args.out+".putative.bacteria.csv","w")
fileOut_blast.write("contig,id,name,identity,alignment_length,contig_length,adjusted_identity\n")


fileOut=open(args.out,"w")
fileOut.write("contig,id,name,number_mismatches,alignment_length,contig_length,adjusted_identity,flag\n")


for i in set(dict_BWA)-set(dict_blast):
    str=','.join(dict_BWA[i])
    id_BWA = dict_BWA[i][1]
    fileOut.write(str+",no-blast")
    fileOut.write("\n")

print "Contigs which are  mapped to custom microbial database and TREE OF LIFE FROM (BLAST)"

both=set(dict_BWA).intersection(set(dict_blast))
print  (both)


blast_bact=set()


#['k81_6', 'JF145555.1', 'Uncultured_bacterium_clone_ncd1640a09c1_16S_ribosomal_RNA_gene_', '99.412', '510.0', '511.0', '99.2174559687']
#['k81_6', 'gb|AF191073|Strain|UNKNOWN-AF191073|Description|Stealth', 'gb_AF191073_Strain_UNKNOWN-AF191073_Description_Stealth_virus_1_clone_3B43_genomic_sequence__Country__ncbiId_AF191073_1_vipr-id_676234', '159', '445', '511', '64.2696629213']

for i in both:
    print "-----"
    print dict_blast[i]
    print dict_BWA[i]	
    id_blast=dict_blast[i][1].split(".")[0]
    print dict_BWA[i]
    #if args.o=="fungi" or args.o=="protozoa":
    id_BWA=dict_BWA[i][1]
    #elif args.o=="virus":
        #id_BWA=dict_BWA[i][1].split("|")[1]
    genus_blast=dict_blast[i][2].split("_")[0]
    genus_BWA = dict_BWA[i][1].split("_")[0]
    identity_blast=float(dict_blast[i][6])
    identity_BWA = float(dict_BWA[i][6])

    print id_blast,id_BWA,genus_blast,genus_BWA,identity_BWA,identity_blast
    #JF145555 gb|AF191073|Strain|UNKNOWN-AF191073|Description|Stealth Uncultured gb 99.2174559687 99.2174559687


    if id_blast==id_BWA: # the same ids
	print "the same ids"
        str = ','.join(dict_BWA[i])
        fileOut.write(str +",confirmed-blast-id")
        fileOut.write("\n")
    elif genus_blast==genus_BWA: # the same genus
	print "# the same genus"
        str = ','.join(dict_BWA[i])
        fileOut.write(str +",confirmed-blast-genus")
        fileOut.write("\n")
    elif identity_BWA>identity_blast:
	print "identity_BWA>identity_blast"
        str = ','.join(dict_BWA[i])
        fileOut.write(str +",better-than-blast")
        fileOut.write("\n")
    else: #BLAST is better => human then filter, or bacteria then report

	

        if genus_blast!="Homo": # if not human => bact or something else. Need more work here
	    print "Bacteria"	
            if compare(dict_blast2[i]):  #make sure this is significntly different from others
                str = ','.join(dict_blast[i])
                fileOut_blast.write(str)
                fileOut_blast.write("\n")






fileOut_blast.close()
fileOut.close()

num_lines_virus = sum(1 for line in open(args.out))
num_lines_bact = sum(1 for line in open(args.out+".putative.bacteria.csv"))

num_lines_virus=num_lines_virus-1
num_lines_bact=num_lines_bact-1


print num_lines_virus, " filtered contigs are saved to", args.out
print num_lines_bact, "contigs which better match bacteria are detected. We report only ones which are different by at least 2%", args.out+".putative.bacteria.csv"


print "done!"




