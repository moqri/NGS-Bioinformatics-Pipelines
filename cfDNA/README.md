https://dnmtools.readthedocs.io/en/latest/

```
module add trim-galore/0.6.7 
```

```
f=SRR15143380
hg38=/labs/mpsnyder/moqri/data_all/ref/hg38/
p=60
trim_galore -q 20 --length 50 --three_prime_clip_R1 20 --clip_R1 20 --clip_R2 20 --paired "$f"_1.fastq.gz "$f"_2.fastq.gz -j $p

```

```
bismark --genome $hg38 --parallel $p -1 "$f"_1_val_1.fq.gz -2 "$f"_2_val_2.fq.gz

/labs/mpsnyder/moqri/soft/Bismark-0.24.1/bismark_methylation_extractor "$f"_1_val_1_bismark_bt2_pe.bam --parallel $p --bedGraph --comprehensive --ucsc

awk '$6>10' "$f"_1_val_1_bismark_bt2_pe.bismark.cov > $f.cov
awk '{print $1,$2,$2+1,$4}' $f.cov > $f.bed

/labs/mpsnyder/moqri/soft/ucsc/bedGraphToBigWig $f.bed ../ch.size $f.bw
cp $f.bw /labs/vsebast/shared/PRC2_Clock/data/cfDNA/
/labs/mpsnyder/moqri/soft/Bismark-0.24.1/bismark2bedGraph CpG_context_SRR15143251_1_val_1_bismark_bt2_pe.t
xt -o bed --ucsc  ample
```

```
bismark_genome_preparation /labs/mpsnyder/moqri/data_all/ref/hg38/ --parallel $p
bismark_methylation_extractor SRR15143251_1_val_1_bismark_bt2_pe.bam --parallel $p --bedGraph --comprehensive --ucsc
dnmtools abismal -i $ind -s "$f".stats "$f"_1_val_1.fq "$f"_2_val_2.fq -t $p -v | samtools view -bS > "$f".bam
```
