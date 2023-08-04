https://dnmtools.readthedocs.io/en/latest/


```
module add trim-galore/0.6.7 
f=SRR15143251
p=60
trim_galore -q 20 --length 50 --three_prime_clip_R1 20 --clip_R1 20 --clip_R2 20 --paired "$f"_1.fastq "$f"_2.fastq -j $p
dnmtools abismal -i $ind -s "$f".stats "$f"_1_val_1.fq "$f"_2_val_2.fq -t $p -v | samtools view -bS > "$f".bam
```


```
bismark_genome_preparation /labs/mpsnyder/moqri/data_all/ref/hg38/ --parallel $p
bismark --genome $hg38_folder --parallel $p -1 SRR15143251_1_val_1.fq -2 SRR15143251_2_val_2.fq
bismark_methylation_extractor SRR15143251_1_val_1_bismark_bt2_pe.bam --parallel $p --bedGraph --comprehensive --ucsc
```
