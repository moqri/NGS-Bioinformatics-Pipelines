GSE186888
```
for f in ./*.beta
do 
 /labs/mpsnyder/moqri/soft/wgbs_tools/wgbstools beta2bed "$f" > "$f".bed
done

for f in ./*.bed
do 
 /labs/mpsnyder/moqri/soft/ucsc/liftOver $f /labs/mpsnyder/moqri/soft/ucsc/hg19ToHg38.over.chain.gz $f.hg38 $f.log
done

for f in hg38/*hg38
do 
	awk '!a[$1 $2]++ { print ;}' $f > "${f:12:3}".bed
done


for f in $(seq 194 194)
do 
 awk '{printf "%s\t%d\t%s\t%s\t%.3f\t%d\n",$1,$2,"+","CpG",$4/$5,$4+$5}' bed/$f.bed > meth/$f.meth
 awk -v OFS="\t" '{print $1, $2, $2+1, $4":"$6, $5, $3}' meth/$f.meth > meth_bed/$f.bed
 cut -f 1-3,5 meth_bed/$f.bed | /labs/mpsnyder/moqri/soft/ucsc/wigToBigWig /dev/stdin /labs/mpsnyder/moqri/soft/ucsc/hg38.chrom.sizes bw/$f.bw
 C_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o meth/$f.sorted meth/sorted/$f.meth
done

scp moqri@smsh11dsu-srcf-d15-35.scg.stanford.edu:/labs/vsebast/shared/PRC2_Clock/data/meth_b/geo/beta/bw/174.bw .

dnmtools merge 172.sorted  174.sorted  176.sorted  178.sorted  180.sorted  182.sorted  184.sorted  186.sorted  188.sorted  190.sorted  192.sorted 173.sorted  175.sorted  177.sorted  179.sorted  181.sorted  183.sorted  185.sorted  187.sorted  189.sorted  191.sorted  193.sorted 194.sorted -o merged.meth -v

dnmtools hmr-rep merged.meth -o merged.hmr
```
