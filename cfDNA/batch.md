```
module spider sratoolkit
cd fastq
for s in SRR15143257 SRR15143260 SRR15143267 SRR15143268 SRR15143291 SRR15143317 SRR15143380 SRR9982183 SRR9982229 SRR9982509 SRR9982839 SRR9982938
do
 mkdir $s
 cd $s
 fasterq-dump $s -v
 cd ..
done
```

```
module add trim-galore/0.6.7 
hg38=/labs/mpsnyder/moqri/data_all/ref/hg38/
cd trim
for f in SRR15143257 SRR15143260 SRR15143267 SRR15143268 SRR15143291 SRR15143317 SRR15143380 SRR9982183 SRR9982229 SRR9982509 SRR9982839 SRR9982938
do
 mkdir $f
 cd $f
 trim_galore -q 20 --length 50 --three_prime_clip_R1 20 --clip_R1 20 --clip_R2 20 --paired ../../fastq/$f/"$f"_1.fastq ../../fastq/$f/"$f"_2.fastq -j $p
 cd ..
done
```

```
module add bowtie2
module add bismark

hg38=/labs/mpsnyder/moqri/data_all/ref/hg38/
p=60
cd map
for f in SRR15143257 SRR15143260 SRR15143267 SRR15143268 SRR15143291 SRR15143317 SRR15143380 SRR9982183 SRR9982229 SRR9982509 SRR9982839 SRR9982938
do
 mkdir $f
 cd $f
 bismark --genome $hg38 --parallel $p -1 ../../trim/$f/"$f"_1_val_1.fq -2 ../../trim/$f/"$f"_2_val_2.fq
 /labs/mpsnyder/moqri/soft/Bismark-0.24.1/bismark_methylation_extractor "$f"_1_val_1_bismark_bt2_pe.bam --parallel $p --bedGraph --comprehensive --ucsc
 gzip -d "$f"_1_val_1_bismark_bt2_pe.bismark.cov.gz
 awk '{print $1,$2,$2+1,$4}' "$f"_1_val_1_bismark_bt2_pe.bismark.cov > $f.bed
 cd ..
done
```



