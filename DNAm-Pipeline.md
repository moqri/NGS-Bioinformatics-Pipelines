```
module add trim-galore/0.6.7 
module load samtools
```
```
f=<fastq filename prefix> # FASTQ for IP or Control
p=<number of processors>
```
```
trim_galore --paired -q 0 --length 0 "$f"_1.fastq.gz  "$f"_2.fastq.gz -j $p
abismal -i $ind "$f"_1_val_1.fq.gz "$f"_2_val_2.fq.gz -t $p -v | samtools view -b > "$f".bam 
format_reads -f abismal "$f".bam -o "$f"_f.sam 
samtools sort -O bam -o "$f"_fs.bam "$f"_f.sam -@ 16 -T tmp # -m 8G ?
duplicate-remover -S "$f"_stat.txt "$f"_fs.bam "$f"_fsd.sam
methcounts -c $genome -o "$f".meth "$f"_fsd.sam -v -n
```
