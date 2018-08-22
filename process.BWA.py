import pysam
import csv
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('bam', help='--')
ap.add_argument('out', help='--')
args = ap.parse_args()



fileOut=open(args.out,"w")


fileOut.write("contig,id,name,number_mismatches,alignment_length,contig_length,adjusted_identity\n")



with pysam.AlignmentFile(args.bam, 'rb', check_sq=False) as input_fo:
            for read in input_fo.fetch():
                number_mismatches = int(read.get_tag('NM'))
                read_length = int(read.infer_read_length())
                alignment_length = int(read.query_alignment_length)
                soft = read_length - alignment_length
                number_mismatches += soft
                
                
                
                adjusted_identity=(1.0-float(number_mismatches/float(alignment_length)))*100
                fileOut.write(read.query_name+","+read.reference_name+","+read.reference_name+","+str(number_mismatches)+","+str(alignment_length)+","+str(read_length)+","+str(adjusted_identity))
                fileOut.write("\n")

fileOut.close()


print "Results are here",args.out
print "Done!"
