#!/bin/bash

############ This section until mapping was done on SCG #################


index=$1
id=`sed "${SLURM_ARRAY_TASK_ID}p;d" ${index} | awk '{print $1}'`

genome_index=$2

workdir=$3

fastas=${workdir}/Fastas


cd $workdir

mkdir Trimmed

cd $fastas/

echo '************************************'
echo Start ${id} Trimming:
echo '************************************'
date

trim_galore --paired --gzip --cores 3 ${id}_*1.fq.gz ${id}_*2.fq.gz --output_dir ../Trimmed/

#SE version
#trim_galore --gzip --cores 3 ${id}.fastq.gz  -o ../Trimmed/


echo '************************************'
echo Start $id Mapping:
echo '************************************'
date

mkdir ../ChIP_Bams

cd ChIP_Bams

R1=${id}_*1*.fq.gz
R2=${id}_*2*.fq.gz

echo Processing $R1 $R2

bowtie2 -p 13 -x $genome_index/GRCh38 -1 $workdir/Trimmed/$R1 -2 $workdir/Trimmed/$R2 -S $workdir/ChIP_Bams/${id}_.sam #genome_index/GRCh38 doesnt refer to GRCh38 as a subfolder, its necessary to id the genome index tho
# SE VERSION bowtie2 -p 15 -x $genome_index/GRCh38  $workdir/Trimmed/$R1 -S $workdir/ChIP_Bams/${id}_.sam #genome_index/GRCh38 doesnt refer to GRCh38 as a subfolder, its necessary to id the genome index tho
samtools view -bS ${id}_.sam > ${id}_.bam #also ran -Su instead of -bS. -Su is to do with converting to a bam, which is what ENCODE does, didn't make a difference however.
samtools sort ${id}_.bam -o ${id}_sorted.bam
samtools flagstat ${id}_sorted.bam > ${id}_stats.txt
samtools index ${id}_sorted.bam ${id}_sorted.bai
# Clean up
#rm -f ${id}_.sam
#rm -f $id.bam




###################### Below is the MACS2 peakcalling, as well as bam filtering, done locally ############################

#getting the t cell data

get /labs/vsebast/shared/PRC2_Clock/data/meth_t/* /Users/dsimps93/Project_PRC2_blood_july2023/PRC2_index/data/t_cell

#Write to an output file I think this has to be ran outside the script. It seems to inhibit the loop from running.
#script ChIP_Fibropass_Err_070323.txt


mkdir filt_BAMs
mkdir macs_out
mkdir macs_sorted
mkdir macs_bed_filt
mkdir macs_bedsort
mkdir Biggywigs_ourMACs

mkdir macs_BEDs

index=`cat samps.txt`

# filtering and calling fragment size to use later in signal calling
for fileID in $index; do

    sambamba view -h -t 15 -f bam -F "[XS] == null and not unmapped and not duplicate and proper_pair and mapping_quality > 35" BAMs/${fileID}_*.bam > filt_BAMs/${fileID}_filt.bam ;
    #take out proper_pair for SE ver

    macs2 predictd -i filt_BAMs/${fileID}_*.bam -g hs -m 5 20

done

#frags sizes in order of samples
frags="297 295 298 296 296 295"



####running CD4s merged

CD4s="EZH2_CD4_d32 EZH2_CD4_d43 EZH2_CD4_d48"

FRAGLEN="297" #mean frag length of three samples

macs2 callpeak -t filt_BAMs/EZH2_CD4_*.bam \
        -c filt_BAMs/Inp1_CD4_*.bam \
        -f BAMPE -g hs -p 1e-2 --keep-dup all -B --SPMR --nomodel --extsize $FRAGLEN \
        -n EZH2_CD4_pooled \
        --outdir macs_out \
    ;

fileID="EZH2_CD4_pooled"

macs2 bdgcmp -t macs_out/${fileID}_treat_pileup.bdg -c macs_out/${fileID}_control_lambda.bdg -o macs_BEDs/${fileID}_fe.bdg -m FE 

#filtering for only autosomes, glitch in MAC2 can create strange chroms that arent present or out of bounds
grep -h -E '^([1-9]|1[0-9]|2[0-3]|X|Y|MT)\b' macs_BEDs/${fileID}_fe.bdg > macs_bed_filt/${fileID}_filt_fe.bdg 
#UCSC ver grep -h -E '^(chr[1-9]|chr1[0-9]|chr2[0-3]|chrX|chrY|chrMT|chrM)\b' macs_BEDs/${fileID}_fe.bdg > macs_bed_filt/${fileID}_filt_fe.bdg 


sort -k1,1 -k2,2n macs_bed_filt/${fileID}_*_fe.bdg > macs_bedsort/${fileID}_sorted_fe.bdg 

bedGraphToBigWig macs_bedsort/${fileID}_*_fe.bdg /Users/dsimps93/genomes/GRCh38/hg38_ensembEd.chrom.sizes.txt Biggywigs_ourMACs/${fileID}_fe.bw 


#adding pvalue- poisson dist here too

#theres a setting I missed in ppois that ENCODE do, gonna include and redo here so that I can see if it makes a difference on the already filtered samp to get the signal we want

# Compute sval = min(no. of reads in ChIP, no. of reads in control) / 1,000,000
chipReads1=`samtools view -c filt_BAMs/${fileID1}.bam`
chipReads2=`samtools view -c filt_BAMs/${fileID2}.bam` 
chipReads=$(awk -v total_reads=$((chipReads1 + chipReads2)) 'BEGIN { printf "%.6f", total_reads / 2 / 1000000 }')

controlReads1=`samtools view -c filt_BAMs/${cnt1}.bam`
controlReads2=`samtools view -c filt_BAMs/${cnt2}.bam` 
controlReads=$(awk -v total_reads=$((controlReads1 + controlReads2)) 'BEGIN { printf "%.6f", total_reads / 2 / 1000000 }')

sval=$(echo "${chipReads} ${controlReads}" | awk '$1>$2{printf "%f",$2} $1<=$2{printf "%f",$1}');


macs2 bdgcmp -t macs_out/${fileID}_treat_pileup.bdg -c macs_out/${fileID}_control_lambda.bdg -o macs_BEDs/${fileID}_sval_ppois.bdg -m ppois -S ${sval}

#grep -h -E '^([1-9]|1[0-9]|2[0-3]|X|Y|MT)\b' macs_BEDs/${fileID}_*ppois.bdg > macs_bed_filt/${fileID}_filt_ppois.bdg 
grep -h -E '^(chr[1-9]|chr1[0-9]|chr2[0-3]|chrX|chrY|chrMT|chrM)\b' macs_BEDs/${fileID}_*sval_ppois.bdg > macs_bed_filt/${fileID}_filt_sval_ppois.bdg 


sort -k1,1 -k2,2n macs_bed_filt/${fileID}_*sval_ppois.bdg > macs_bedsort/${fileID}_sorted_sval_ppois.bdg 

bedGraphToBigWig macs_bedsort/${fileID}_*sval_ppois.bdg /Users/dsimps93/genomes/GRCh38/hg38.UCSCchrom.sizes Biggywigs_ourMACs/${fileID}_sval_ppois.bw 



