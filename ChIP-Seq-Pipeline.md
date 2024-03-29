## Mapping
Map IP and Control samples

```
module add sratoolkit
module add trim-galore/0.6.7 
module add bowtie2 
module add samtools
module add sambamba
module add macs2   

f=<fastq filename prefix> # FASTQ for IP or Control
p=<number of processors>
index=<bowtie2 index path>

fasterq-dump $f -p "-e$p"
trim_galore --paired -q 0 --length 0 "$f"_R1_001.fastq.gz  "$f"_R2_001.fastq.gz -j $p
bowtie2 -q -x $index -1 "$f"_R1_001_val_1.fq.gz -2 "$f"_R1_001_val_1.fq.gz -S $f.sam --local --no-unal --very-sensitive -X 2000 -p $p
### faster less sensitive version: 
bowtie2 -q -x $index -1 "$f"_R1_001_val_1.fq.gz -2 "$f"_R1_001_val_1.fq.gz -S $f.sam -p $p

samtools view -h -S -b -o $f.bam $f.sam -@ $p
sambamba sort -o "$f"_s.bam $f.bam -t $p
sambamba view -h -t $p -f bam -F "[XS] == null and not unmapped  and not duplicate" "$f"_s.bam > "$f"_d.bam

```

## Peak-calling
```
c=<control prefix>
macs2 callpeak -t "$f"_d.bam -c "$c"_d.bam -f BAMPE -g hs -B -n $f --outdir $f
bedGraphToBigWig $f.bdg chrom.sizes $f.bw
```
