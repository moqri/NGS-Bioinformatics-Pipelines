## Mapping
```
module add trim-galore/0.6.7 
module add bowtie2 
module load samtools
module load sambamba

f=<fastq filename prefix> # FASTQ for IP or Control
p=<number of processors>
index=<bowtie2 index path>

fasterq-dump $f -e60
trim_galore --paired -q 0 --length 0 "$f"_R1_001.fastq.gz  "$f"_R2_001.fastq.gz -j $p
bowtie2 -q -x $index -1 "$f"_R1_001_val_1.fq.gz -2 "$f"_R1_001_val_1.fq.gz -S $f.sam --local --no-unal --very-sensitive -X 2000 -p $p

samtools view -h -S -b -o $f.bam $f.sam -@ $p
sambamba sort -o "$f"_s.bam $f.bam -t $p
```

## Peack calling
```
Genrich  -t "$s"_n.bam -o $s.bed -j  -y  -r  -e chrM  -v
```
## Consencus
```
#concensus
```
