ChIP: 
https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq

```
module add samtools
module add picard/3.0.0
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
