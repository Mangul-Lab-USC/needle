import pysam
import csv
import argparse
import os
import gzip

ap = argparse.ArgumentParser()
ap.add_argument('bam', help='--')
ap.add_argument("-o",help="[virus, fungi or protozoa]. By default virus", default='virus',type=str)
ap.add_argument('out', help='--')
args = ap.parse_args()



dict={}
DIR=os.path.dirname(os.path.realpath(__file__))

if args.o=="virus":

    #>gb|AB261990|Strain|M2|Description|Lymphocytic choriomeningitis virus genomic RNA, segment S, nearly complete sequence including cds, strain: M2.|Country||ncbiId|AB261990.1|vipr-id|284567
    print "Open look up table VIPR.table.txt.gz"
    
    with gzip.open(DIR+'/VIPR.table.txt.gz', 'rb') as f:
        for line in f:
            id=line.split("|")[1]
            name=line.split("|")[4].replace(">","").replace(",","").replace(" ","_").replace("|","_").replace(".","_").replace(":","").replace("\n","")
            dict[id]=name
elif args.o=="protozoa":
    #>NC_001905.3_Leishmania_major
    #BAM : NC_001905.3
    print "Open look up table protozoa.table.txt"
    file=open(DIR+'/protozoa.table.csv')
    reader=csv.reader(file)
    for line in reader:
            id=line[0]
            name='_'.join(line[1:])
            dict[id]=name
    file.close()

elif args.o=="fungi":
    #>NC_001905.3_Leishmania_major
    #BAM : NC_001905.3
    print "Open look up table  fungi.table.txt"
    file=open(DIR+'/fungi.table.csv')
    reader=csv.reader(file)
    for line in reader:
        id=line[0]
        name='_'.join(line[1:])
        dict[id]=name
    file.close()
else:
    print "Error. Needs to be : virus, fungi, protozoa"
    sys.exit(0)






fileOut=open(args.out,"w")


fileOut.write("contig,id,name,number_mismatches,alignment_length,contig_length,adjusted_identity\n")




with pysam.AlignmentFile(args.bam, 'rb', check_sq=False) as input_fo:
            for read in input_fo.fetch(until_eof=True):
                number_mismatches = int(read.get_tag('NM'))
                read_length = int(read.infer_read_length())
                alignment_length = int(read.query_alignment_length)
                soft = read_length - alignment_length
                number_mismatches += soft
                
                if args.o=="virus":
                    id_BWA=read.reference_name.split("|")[1]
                elif args.o=="protozoa" or args.o=="fungi":
                    id_BWA=read.reference_name
                
                
                adjusted_identity=(1.0-float(number_mismatches/float(alignment_length)))*100
                fileOut.write(read.query_name+","+id_BWA+","+dict[id_BWA]+","+str(number_mismatches)+","+str(alignment_length)+","+str(read_length)+","+str(adjusted_identity))
                fileOut.write("\n")

fileOut.close()


print "Results are here",args.out
print "Done!"
