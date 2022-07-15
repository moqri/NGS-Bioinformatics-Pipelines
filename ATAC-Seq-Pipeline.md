Source: 
* https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/
* https://informatics.fas.harvard.edu/atac-seq-guidelines.html

## Set up
* You will need [Genrich](https://github.com/jsh58/Genrich) (if using Genrich for peack-calling) and [MSPS](https://genometric.github.io/MSPC/) if mergeing peacks.

```
module add sratoolkit
module add trim-galore/0.6.7 
module add bowtie2 
module add samtools
module add sambamba
module add macs3
```
## Trimming
```
f=<fastq filename prefix> # e.g. SRA ID
p=64 #number of processors
ind=/mm10/mm10 # bowtie2 index path to .bt files

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
bowtie2 -q -x $ind -1 "$f"_1_val_1.fq -2 "$f"_2_val_2.fq -S $f.sam --local --no-unal --very-sensitive -X 2000 -p $p
samtools view -bS $f.sam > $f.bam -@ $p
samtools view -h -q 30 $f.bam -b > $f.q.bam -@ $p
```
## Peack calling using MACS3
### Single reads
```
samtools sort $f.q.bam > "$f".s.bam -@ 16 
java -Xmx120g -jar picard.jar MarkDuplicates I=$f.s.bam O=$f.pic.bam M=tmp REMOVE_DUPLICATES=true
macs3 callpeak -t $f.pic.bam -g hs -n $f.bed -q 0.01 --broad
```
### Paired reads
```
samtools view -h -b -F 1804 -f 2 $f.q.bam -b > $f.qq.bam -@ $p 
samtools view -hf 0x2 $f.qq.bam -b > $f.p.bam #proper paired
samtools sort $f.p.bam > "$f".s.bam -@ 16
java -Xmx120g -jar picard.jar MarkDuplicates I=$f.s.bam O=$f.pic.bam M=tmp REMOVE_DUPLICATES=true
macs3 callpeak -f BAMPE -t $f.pic.bam -g hs -n $f.bed -q 0.01 --broad
```
## Peack calling using Genrich
```
samtools sort -n $f.bam > "$f"_n.bam -@ 16
Genrich  -t "$f"_n.bam -o $f.bed -j  -y  -r  -e chrM  -v
```
## Consencus
```
mspc -i *.bed -r bio -w 1e-4 -s 1e-8
```
