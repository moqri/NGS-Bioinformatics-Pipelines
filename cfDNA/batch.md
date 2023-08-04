```
module spider sratoolkit
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

for f in SRR15143257 SRR15143260 SRR15143267 SRR15143268 SRR15143291 SRR15143317 SRR15143380 SRR9982183 SRR9982229 SRR9982509 SRR9982839 SRR9982938
do
 hg38=/labs/mpsnyder/moqri/data_all/ref/hg38/
 mkdir $f
 cd $f
 trim_galore -q 20 --length 50 --three_prime_clip_R1 20 --clip_R1 20 --clip_R2 20 --paired ../../fastq/$f/"$f"_1.fastq ../../fastq/$f/"$f"_2.fastq -j $p
 cd ..
done
```

