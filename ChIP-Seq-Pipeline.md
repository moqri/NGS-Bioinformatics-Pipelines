
```
module add trim-galore/0.6.7 
module add bowtie2 
module load samtools
module load sambamba

f=<fastq filename>
p=<number of processors>
index=<bowtie2 index path>

trim_galore --paired -q 0 --length 0 "$f"_R1_001.fastq.gz  "$f"_R2_001.fastq.gz -j $p
bowtie2 -q -x $index -1 "$f"_R1_001_val_1.fq.gz -2 "$f"_R1_001_val_1.fq.gz -S $f.sam --local --no-unal --very-sensitive -X 2000 -p $p

samtools view -h -S -b -o $f.bam $f.sam -@ $p
sambamba sort -o "$f"_s.bam $f.bam -t $p
sambamba view -h -t $p -f bam -F "[XS] == null and not unmapped  and not duplicate" "$f"_s.bam > "$f"_d.bam

```

Do the same for the input $c

```
module add macs2   
macs2 callpeak -t "$f"_d.bam -c "$c"_f.bam -f BAMPE -g hs -B -n $f --outdir macs2
```
