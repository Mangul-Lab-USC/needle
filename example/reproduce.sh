../tools/MiniConda/bin/python ../process.BWA.py toy.example.virus.megahit.contigs.SV.bam toy.example.virus.megahit.contigs.SV.filtered.csv
python ../tree.of.life.filter.py ETAG001.blast.csv Sample_ETAG001.Aligned.out_after_rRNA.fasta.virus.megahit.contigs.csv ETAG001.final.csv

python ../tree.of.life.filter.py toy.example.cat.unmapped.fastq.fungi.megahit.contigs.BLAST.house.format.csv toy.example.cat.unmapped.fastq.fungi.megahit.contigs.SV.csv toy.example.cat.unmapped.fungi.filtered.csv

python ../tree.of.life.filter.py toy.example.virus.megahit.contigs.BLAST.house.format.csv toy.example.virus.megahit.contigs.SV.csv toy.example.virus.megahit.contigs.SV.filtered.csv
