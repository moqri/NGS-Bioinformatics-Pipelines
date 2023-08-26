GSE186888
```
for f in beta/*.beta
do 
 /labs/mpsnyder/moqri/soft/wgbs_tools/wgbstools beta2bed "$f" > bed/$f.bed
done

for f in *.bed
do 
 /labs/mpsnyder/moqri/soft/ucsc/liftOver $f /labs/mpsnyder/moqri/soft/ucsc/hg19ToHg38.over.chain.gz ../hg38/$f.hg38 ../log/$f.log
 awk '!a[$1 $2]++ { print ;}' ../hg38/$f.hg38 > ../"${f:7:3}".bed
done

for f in $(seq 172 194)
do 
 awk '{printf "%s\t%d\t%s\t%s\t%.3f\t%d\n",$1,$2,"+","CpG",$4/$5,$4+$5}' bed/$f.bed > meth/$f.meth
 awk -v OFS="\t" '{print $1, $2, $2+1, $4":"$6, $5, $3}' meth/$f.meth > meth_bed/$f.bed
 cut -f 1-3,5 meth_bed/$f.bed | /labs/mpsnyder/moqri/soft/ucsc/wigToBigWig /dev/stdin /labs/mpsnyder/moqri/soft/ucsc/hg38.chrom.sizes bw/$f.bw
 C_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o meth/sorted/$f.sorted meth/$f.meth
done

dnmtools merge 172.sorted  174.sorted  176.sorted  178.sorted  180.sorted  182.sorted  184.sorted  186.sorted  188.sorted  190.sorted  192.sorted 173.sorted  175.sorted  177.sorted  179.sorted  181.sorted  183.sorted  185.sorted  187.sorted  189.sorted  191.sorted  193.sorted 194.sorted -o merged.meth -v

dnmtools hmr-rep merged.meth -o merged.hmr
```

ELSA-SEQ
```
dnmtools abismal -i /labs/vsebast/moqri/epi10x/data/neo_hmr.idx -o map/$f.sam trim/$f/"$f"_1_val_1.fq trim/$f/"$f"_2_val_2.fq -t $t

samtools sort -O bam -o map/$f.bam map/$f.sam -@ 16 -T tmp 

dnmtools counts -c /labs/vsebast/moqri/epi10x/data/neo_hmr.fasta map/$f.bam -n -o meth/$f.meth -v

dnmtools sym meth/$f.meth -o meth/$f.sym

awk '{split($1,a,":"); print a[1],a[2]+$2,$3,$4,$5,$6}' meth/$f.sym > meth/$f.hg38
awk '{print $1,$2,$2+1,$5}' meth/$f.hg38 | /labs/mpsnyder/moqri/soft/ucsc/wigToBigWig /dev/stdin /labs/mpsnyder/moqri/soft/ucsc/hg38.chrom.sizes bw/$f.bw
```

GSE117593 Willson
```
module add trim-galore/0.6.7 
p=64
t=64

for i in $(seq 22 30)
do 
 f="SRR75870"$i
 trim_galore -q 0 --length 0 fastq/"$f".fastq.gz -o trim -j $p
 dnmtools abismal -i /labs/vsebast/moqri/epi10x/data/neo_hmr.idx -o map/$f.sam trim/"$f"_trimmed.fq.gz -t $t
 samtools sort -O bam -o map/$f.bam map/$f.sam -@ 16 -T tmp 
 dnmtools counts -c /labs/vsebast/moqri/epi10x/data/neo_hmr.fasta map/$f.bam -n -o meth/$f.meth -v
 dnmtools sym meth/$f.meth -o meth/$f.sym
 awk '{split($1,a,":"); print a[1],a[2]+$2,$3,$4,$5,$6}' meth/$f.sym > meth/$f.hg38 awk '{print $1,$2,$2+1,$5}' meth/$f.hg38 | /labs/mpsnyder/moqri/soft/ucsc/wigToBigWig /dev/stdin /labs/mpsnyder/moqri/soft/ucsc/hg38.chrom.sizes bw/$f.bw
done
```
