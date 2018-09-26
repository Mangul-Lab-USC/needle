#while read line;do mv ../${line}* ./;done</u/home/s/serghei/code/seeing.beyond.target.results/outcomes/WXS.fileNames.txt

ls *bam | awk -F ".bam" '{print $1}' >samples.txt

while read line
do
echo "/u/home/s/serghei/project/code/needle/needle.sh ${line}.bam ${line}">run.${line}.sh
qsub -cwd -V -N  needle.${line} -l h_data=12G,highp,time=10:00:00 run.${line}.sh
done<samples.txt 
