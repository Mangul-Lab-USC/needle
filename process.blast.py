import argparse
import re
import sys
import csv


ap = argparse.ArgumentParser()
ap.add_argument('input_short', help='--')
ap.add_argument('input_long', help='--')
ap.add_argument('out', help='---')
args = ap.parse_args()




file = open(args.input_short,"r")
file_long=open(args.input_long,"r")
fileOut=open(args.out, "w")




#open long
#AB261990.1  Lymphocytic choriomeningitis virus genomic RNA, segme...  1692    0.0

dict={}

reader=csv.reader(file_long,delimiter=' ')
for line in reader:
    if len(line)>0:
        if len(line[0])>0:
            if line[0][0]=='>':
                id=line[0].split('>')[1]
                name='_'.join(line[1:]).replace(",","")
                dict[id]=name
file_long.close()




fileOut.write("contig,id,name,identity,alignment_length,contig_length,adjusted_identity\n")

#k21_2    gi|108743534|dbj|AB261990.1|    100.000    916    916    0

for line in file:
    
    if re.match("# Fields: query id, subject id, % identity, query length, alignment length, mismatches",line):
        new_column = ', adjusted_identity'
        fields = line.split()
        del fields[0:2]
        new_fields = " "
        new_fields = new_fields.join(fields)
        new_fields = new_fields.rstrip() + new_column
    
    elif(line[0]!='#'):
            split= line.split()
            query_length = float(split[3])
            identity = float(split[2])
            alignment_length = float(split[4])
            adjusted_identity = str(alignment_length * identity / query_length)
            contig=split[0]
            id=split[1].split('|')[3]
            fileOut.write(contig+","+id+","+dict[id]+","+str(identity)+","+str(alignment_length)+","+str(query_length)+","+str(adjusted_identity))
            fileOut.write("\n")
            
            
            new_identity = ","
            new_identity = new_identity.join(split)
            new_identity = new_identity.rstrip() + ","+ adjusted_identity

file.close()













