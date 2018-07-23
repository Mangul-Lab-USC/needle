#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('bam')
parser.add_argument('outdir')
parser.add_argument('-f', '--force', action='store_true', default=False,
                    help='Forse [default %(default)s]')
EOF

DIR_CODE=`dirname $(readlink -f "$0")`

echo required infile: "$INBAM"
echo required outfile: "$OUTDIR"


ORGANISM='human'
DB="$DIR_CODE/db_$ORGANISM"

#Add MiniConda to PATH if it's available.
if [ -d "$DIR_CODE/tools/MiniConda/bin" ]; then
    echo "Add MiniConda to PATH if it's available"
    export PATH="$DIR_CODE/tools/MiniConda/bin:$PATH"
fi

#Convert to absolute paths.
BAM=`readlink -m "$BAM"`
OUTDIR=`readlink -m "$OUTDIR"`

#Check if BAM exists.
if [ ! -e "$BAM" ]
then
    echo "Error: $BAM doesn't exist." >&2
    exit 1
fi

#Check if OUTDIR exists, then make it.
if [ -d "$OUTDIR" ]
then
    if [[ $FORCE ]]
    then
        rm -fr "$OUTDIR"
    else
        echo "Error: The directory $OUTDIR exists. Please choose a" \
            'different directory in which to save results of the analysis, or' \
            'use the -f option to overwrite the directory.' >&2
        exit 1
    fi
fi
mkdir -p "$OUTDIR"

start=`date +%s`
echo  "Start needle analysis ... "$start

megahit=${DIR_CODE}/tools/megahit/megahit
prefix=$(basename $BAM | awk -F ".bam" '{print $1}')
SAMPLE=${OUTDIR}"/"${prefix}



echo "Extract unmapped reads from " $BAM
samtools view -f 0x4 -bh $BAM | samtools bam2fq - >${SAMPLE}.unmapped.fastq
#samtools view -bh $BAM NC_007605 | samtools fastq - > ${SAMPLE}.NC_007605.fastq
#rm -fr ${SAMPLE}.NC_007605.fastq
#cat ${SAMPLE}.unmapped.fastq ${SAMPLE}.NC_007605.fastq>${SAMPLE}.cat.unmapped.fastq
#rm -fr ${SAMPLE}.unmapped.fastq
UNMAPPED=${SAMPLE}.unmapped.fastq



bwa mem -a ${DB}/viral.vipr/NONFLU_All.fastq $UNMAPPED | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.virus.bam
bwa mem -a ${DB}/fungi/fungi.ncbi.february.3.2018.fasta $UNMAPPED | samtools view -S -b -F 4 - |  samtools sort - >${SAMPLE}.fungi.bam
bwa mem -a ${DB}/protozoa/protozoa.ncbi.february.3.2018.fasta $UNMAPPED | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.protozoa.bam

samtools index ${SAMPLE}.virus.bam
samtools index ${SAMPLE}.fungi.bam
samtools index ${SAMPLE}.protozoa.bam
samtools fastq ${SAMPLE}.virus.bam >${SAMPLE}.virus.fastq
samtools fastq ${SAMPLE}.fungi.bam >${SAMPLE}.fungi.fastq
samtools fastq ${SAMPLE}.protozoa.bam >${SAMPLE}.protozoa.fastq

rm -fr ${SAMPLE}*bam
rm -fr ${SAMPLE}*bai

$megahit --k-step 10 -r ${SAMPLE}.virus.fastq -o ${SAMPLE}.virus.megahit --out-prefix virus.megahit
$megahit --k-step 10 -r ${SAMPLE}.fungi.fastq -o ${SAMPLE}.fungi.megahit --out-prefix fungi.megahit
$megahit --k-step 10 -r ${SAMPLE}.protozoa.fastq -o ${SAMPLE}.protozoa.megahit --out-prefix protozoa.megahit
mv ${SAMPLE}.virus.megahit/virus.megahit.contigs.fa ${SAMPLE}.virus.megahit.contigs.fa
mv ${SAMPLE}.fungi.megahit/fungi.megahit.contigs.fa ${SAMPLE}.fungi.megahit.contigs.fa
mv ${SAMPLE}.protozoa.megahit/protozoa.megahit.contigs.fa ${SAMPLE}.protozoa.megahit.contigs.fa



echo "Map assembled contigs onto human reference"


bwa mem $DB/BWA.index/genome.fa ${SAMPLE}.virus.megahit.contigs.fa | samtools view -bS -F 4 -  >${SAMPLE}.virus.megahit.contigs.human.bam
bwa mem $DB/BWA.index/genome.fa ${SAMPLE}.fungi.megahit.contigs.fa | samtools view -bS -F 4 -  >${SAMPLE}.fungi.megahit.contigs.human.bam
bwa mem $DB/BWA.index/genome.fa ${SAMPLE}.protozoa.megahit.contigs.fa | samtools view -bS -F 4 ->${SAMPLE}.protozoa.megahit.contigs.human.bam


samtools view ${SAMPLE}.virus.megahit.contigs.human.bam | awk '{print $1}'  | sort | uniq >${SAMPLE}.virus.megahit.contigs.human.txt
samtools view ${SAMPLE}.fungi.megahit.contigs.human.bam | awk '{print $1}'  | sort | uniq >${SAMPLE}.fungi.megahit.contigs.human.txt
samtools view ${SAMPLE}.protozoa.megahit.contigs.human.bam | awk '{print $1}'  | sort | uniq >${SAMPLE}.protozoa.megahit.contigs.human.txt


python ${DIR_CODE}/filter.contigs.py ${SAMPLE}.virus.megahit.contigs.fa ${SAMPLE}.virus.megahit.contigs.human.txt ${SAMPLE}.virus.megahit.contigs.filtered.fa
python ${DIR_CODE}/filter.contigs.py ${SAMPLE}.fungi.megahit.contigs.fa ${SAMPLE}.fungi.megahit.contigs.human.txt ${SAMPLE}.fungi.megahit.contigs.filtered.fa
python ${DIR_CODE}/filter.contigs.py ${SAMPLE}.protozoa.megahit.contigs.fa ${SAMPLE}.protozoa.megahit.contigs.human.txt ${SAMPLE}.protozoa.megahit.contigs.filtered.fa

# remove original contigs
rm -fr ${SAMPLE}.virus.megahit.contigs.fa
rm -fr ${SAMPLE}.fungi.megahit.contigs.fa
rm -fr ${SAMPLE}.protozoa.megahit.contigs.fa

#rename filtered contigs file to be as the original ones
mv ${SAMPLE}.virus.megahit.contigs.filtered.fa ${SAMPLE}.virus.megahit.contigs.fa
mv ${SAMPLE}.fungi.megahit.contigs.filtered.fa ${SAMPLE}.fungi.megahit.contigs.fa
mv ${SAMPLE}.protozoa.megahit.contigs.filtered.fa ${SAMPLE}.protozoa.megahit.contigs.fa


# index contigs and map reads onto contigs
bwa index ${SAMPLE}.virus.megahit.contigs.fa
bwa mem ${SAMPLE}.virus.megahit.contigs.fa ${SAMPLE}.virus.fastq | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.megahit.contigs.virus.bam



samtools depth ${SAMPLE}.megahit.contigs.virus.bam>${SAMPLE}.megahit.contigs.virus.cov
samtools view -H ${SAMPLE}.megahit.contigs.virus.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.megahit.contigs.virus.bam | grep -v -e 'XA:Z:' -e 'SA:Z:'| cat ${OUTDIR}/header.sam - | samtools view -b - | samtools depth - >${SAMPLE}.megahit.contigs.virus.uniq.cov




#fungi----
bwa index ${SAMPLE}.fungi.megahit.contigs.fa
bwa mem ${SAMPLE}.fungi.megahit.contigs.fa ${SAMPLE}.fungi.fastq | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.megahit.contigs.fungi.bam


samtools depth ${SAMPLE}.megahit.contigs.fungi.bam>${SAMPLE}.megahit.contigs.fungi.cov
samtools view -H ${SAMPLE}.megahit.contigs.fungi.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.megahit.contigs.fungi.bam | grep -v -e 'XA:Z:' -e 'SA:Z:'| cat ${OUTDIR}/header.sam - | samtools view -b - | samtools depth - >${SAMPLE}.megahit.contigs.fungi.uniq.cov


#protozoa----
bwa index ${SAMPLE}.protozoa.megahit.contigs.fa
bwa mem ${SAMPLE}.protozoa.megahit.contigs.fa ${SAMPLE}.protozoa.fastq | samtools view -S -b -F 4 - | samtools sort - >${SAMPLE}.megahit.contigs.protozoa.bam


samtools depth ${SAMPLE}.megahit.contigs.protozoa.bam>${SAMPLE}.megahit.contigs.protozoa.cov
samtools view -H ${SAMPLE}.megahit.contigs.protozoa.bam >${OUTDIR}/header.sam
samtools view -F 4  ${SAMPLE}.megahit.contigs.protozoa.bam | grep -v -e 'XA:Z:' -e 'SA:Z:'| cat ${OUTDIR}/header.sam - | samtools view -b - | samtools depth - >${SAMPLE}.megahit.contigs.protozoa.uniq.cov

echo "-----------------------------------------------------"
echo "Map assembled contigs onto the microbial references"

bwa mem -a ${DB}/viral.vipr/NONFLU_All.fastq ${SAMPLE}.virus.megahit.contigs.fa | samtools view -bS -F 4 -  >${SAMPLE}.virus.megahit.contigs.SV.bam
bwa mem -a ${DB}/fungi/fungi.ncbi.february.3.2018.fasta ${SAMPLE}.fungi.megahit.contigs.fa  | samtools view -bS -F 4 -  >${SAMPLE}.fungi.megahit.contigs.SV.bam
bwa mem -a ${DB}/protozoa/protozoa.ncbi.february.3.2018.fasta ${SAMPLE}.protozoa.megahit.contigs.fa  | samtools view -bS -F 4 ->${SAMPLE}.protozoa.megahit.contigs.SV.bam





echo "Success!!!"




