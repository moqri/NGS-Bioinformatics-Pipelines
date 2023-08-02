f=SRR15143251
p=64
trim_galore -q 20 --length 50 --three_prime_clip_R1 20 --clip_R1 20 --clip_R2 20 --paired "$f"_1.fastq "$f"_2.fastq -j $p
dnmtools abismal -i $ind -s "$f".stats "$f"_1_val_1.fq "$f"_2_val_2.fq -t $p -v | samtools view -bS > "$f".bam
