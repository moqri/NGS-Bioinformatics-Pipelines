```
f=<fastq filename>
p=<number of processors>
index=<bowtie2 index path>
module add trim-galore/0.6.7
module spider bowtie2
trim_galore --paired -q 0 --length 0 "$f"_R1_001.fastq.gz  "$f"_R2_001.fastq.gz -j $p
bowtie2 -q -x $index -1 "$f"_R1_001_val_1.fq.gz -2 "$f"_R1_001_val_1.fq.gz -S $f.sam --local --no-unal --very-sensitive -X 2000 -p $p
```
