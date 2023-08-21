ChIP: 
https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq

```
module add trim_galore
module add bowtie2
module add samtools
module add picard/3.0.0
```

```
trim_galore --paired fastq/"$f"_1.fq.gz fastq/"$f"_2.fq.gz -o trimmed -j 60
bwt2_idx=/labs/mpsnyder/moqri/data_all/ref/bw2/hg38/GRCh38_noalt_as
bowtie2 -X2000 -x $bwt2_idx -1 trimmed/"$f"_1_val_1.fq.gz -2 trimmed/"f"_2_val_2.fq.gz 2> $f.log -p 16
```

```
samtools view -F 1804 -f 2 -q 30 -u $f.sam -o $f.filtered -O SAM -@16
samtools sort -n $f.filtered -o $f.sorted -O SAM -@16
samtools fixmate -r $f.sorted $f.fixmate -O SAM -@16
samtools view -F 1804 -f 2 -u $f.fixmate -o $f.fixmate.filtered -O SAM -@16 
samtools sort $f.fixmate.filtered -o $f.filtered.sorted -O SAM -@16

picard MarkDuplicates -I $f.filtered.sorted -O $f.marked -METRICS_FILE $f.qc -VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false

samtools view -F 1804 -f 2 $f.marked -o $f.dedup -O SAM -@16
```
