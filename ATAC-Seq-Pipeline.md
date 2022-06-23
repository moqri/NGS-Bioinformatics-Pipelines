Source: 
* https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/
* https://informatics.fas.harvard.edu/atac-seq-guidelines.html

## Set up
* You will need [Genrich](https://github.com/jsh58/Genrich) for peack-calling and [MSPS](https://genometric.github.io/MSPC/) if mergeing peacks.

```
module add sratoolkit
module add trim-galore/0.6.7 
module add bowtie2 
module add samtools
module add sambamba
```
## Trimming
```
f=<fastq filename prefix> # FASTQ for IP or Control
p=<number of processors>
index=<bowtie2 index path>

fasterq-dump $f -p "-e$p"
```
## Mapping
### Single reads
```
trim_galore -q 0 --length 0 "$f".fastq -j $p
bowtie2 -q -x $index -U "$f"_trimmed.fq  -S $f.sam --local --no-unal --very-sensitive -p $p
```
### Paied reads
```
trim_galore --paired -q 0 --length 0 "$f"_1.fastq  "$f"_2.fastq -j $p
bowtie2 -q -x $index -1 "$f"_1_val_1.fq -2 "$f"_2_val_2.fq -S $f.sam --local --no-unal --very-sensitive -X 2000 -p $p
```
## Peack calling
```
samtools view -bS $f.sam > $f.bam -@ $p
samtools sort -n $f.bam > "$f"_n.bam -@ 16
Genrich  -t "$f"_n.bam -o $f.bed -j  -y  -r  -e chrM  -v
```
## Consencus
```
mspc -i *.bed -r bio -w 1e-4 -s 1e-8
```
