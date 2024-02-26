## Setup
```
module add sratoolkit
module add trim-galore/0.6.7 
module add gcc
module add htslib
module add samtools
```
```
f=<fastq filename prefix>
p=60 #number of processors
ind=hg38.abismalidx #abismal index file
ref_genome=hg38.fa #reference genome
```
## Mapping
### from SRA
```
fasterq-dump $f -p "-e$p"
trim_galore --paired -q 0 --length 0 "$f"_1.fastq.gz  "$f"_2.fastq.gz -j $p
abismal -i $ind "$f"_1_val_1.fq.gz "$f"_2_val_2.fq.gz -t $p -v | samtools view -b > "$f".bam 
```
### from sequencing companies
```
trim_galore --paired -q 0 --length 0 "$f"_R1_001.fastq.gz  "$f"_R2_001.fastq.gz -j $p
abismal -i $ind "$f"_R1_001_val_1.fq.gz "$f"_R2_001_val_2.fq.gz -t $p -v | samtools view -b > "$f".bam 
```

## Calculating DNAm levels
```
format_reads -f abismal "$f".bam -o "$f"_f.sam 
samtools sort -O bam -o "$f"_fs.bam "$f"_f.sam -@ 16 -T tmp # -m 8G 
duplicate-remover -S "$f"_stat.txt "$f"_fs.bam "$f"_fsd.sam
methcounts -c $ref_genome -o "$f".meth "$f"_fsd.sam -v -n
symmetric-cpgs -o "$f"_s.meth $f.meth
```
From Bismark BedGraph:
grep -v '_' BedGraph | grep 'chr[1-9][0-9]\?' # only somatic Chrs

