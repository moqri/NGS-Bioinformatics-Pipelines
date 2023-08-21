ENCODE ChIP Doc: 
https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq

```
module add trim_galore
module add bowtie2
module add samtools
module add picard/3.0.0
module add macs2
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

picard MarkDuplicates -I $f.filtered.sorted -O $f.marked.sam -METRICS_FILE $f.qc -VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false

samtools view -F 1804 -f 2 $f.marked.sam -o $f.dedup -O SAM -@16

macs2 predictd -i IP_N2P1/IP_N2P1.dedup
#extsize=297
macs2 callpeak -t IP_N2P1/IP_N2P1.dedup -c Inp_N2P1/Inp_N2P1.dedup  -B --nomodel --extsize 297 --SPMR -n N2P1
macs2 bdgcmp -t N2P1_treat_pileup.bdg -c N2P1_control_lambda.bdg -o N2P1.bdg -m logLR -p 0.00001

LC_COLLATE=C
sort -k1,1 -k2,2n N2P1.bdg > N2P1_s.bdg
awk '$1 !~ /_/ && $1 !~ /M/' N2P1_s.bdg > N2P1_sf.bdg
./bedGraphToBigWig N2P1_sf.bdg hg38.chrom.sizes N2P1.bw
```

```
scp moqri@smsh11dsu-srcf-d15-35.scg.stanford.edu:/oak/stanford/scg/lab_vsebast/shared/fibroblast/vsebast/N2P1.bw .
```
