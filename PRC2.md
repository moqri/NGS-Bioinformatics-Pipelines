ChIP: 
https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq

```
trim_galore --paired fastq/"$f"_1.fq.gz fastq/"$f"_2.fq.gz -o trimmed -j 60
bwt2_idx=/labs/mpsnyder/moqri/data_all/ref/bw2/hg38/GRCh38_noalt_as
bowtie2 -X2000 -x $bwt2_idx -1 trimmed/"$f"_1_val_1.fq.gz -2 trimmed/"f"_2_val_2.fq.gz 2> $f.log -p 16
```
