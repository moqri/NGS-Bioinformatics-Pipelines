```
f=<fastq_filename>
p=<number of processors>
module add trim-galore/0.6.7  
trim_galore --paired -q 0 --length 0 "$f"_R1_001.fastq.gz  "$f"_R2_001.fastq.gz -j $p
```
