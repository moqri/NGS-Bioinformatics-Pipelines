#!/bin/bash

#SBATCH --job-name=AbisIndex
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=250GB
#SBATCH --partition=batch
#SBATCH --account=vsebast
#SBATCH --time=48:00:00
#SBATCH --array=1-6%6

#1-6%2 normally

#for mapping, originally set to 25 cores, 500 gb mem, running two at the same time. Need to rerun without mapping so gonna reduce the resourecs a bit

# sbatch Fibs_DNMTools_Map.sh ListOfSamps.txt /labs/vsebast/DJS/Project_Fibro_Passages/WGBS_DNMTools

module load dnmtools/1.2.2
module load bowtie2/2.5.0
module load samtools/1.9
module load ucsc_tools/309

index=$1
fileID=`sed "${SLURM_ARRAY_TASK_ID}p;d" ${index} | awk '{print $1}'`

bisIndex="/labs/vsebast/DJS/genomes/HomoSap/MethPipe/UCSC_Index/UCSC_hg38.abismalidx"

Index="/labs/vsebast/DJS/genomes/HomoSap/UCSC_ver/hg38.fa"

workdir=$2



echo '************************************'
echo Start ${fileID} Abismal alignment:
echo '************************************'
date

cd $workdir

mkdir  $workdir/Abis_Mapping

# run alignment, adding sorting and conversion to bam in one step
dnmtools abismal -i $bisIndex ../WGBS/Trimmed/${fileID}*1.fq.gz ../WGBS/Trimmed/${fileID}*2.fq.gz -t 25 -v | samtools sort -O bam -@ 25 -o $workdir/Abis_Mapping/${fileID}_abis_mapped_sorted.bam


echo '************************************'
echo Start ${fileID} Uniqing:
echo '************************************'
date


mkdir $workdir/uniqueBams


#Uniqing
dnmtools uniq -v $workdir/Abis_Mapping/${fileID}_abis_mapped_sorted.bam $workdir/uniqueBams/${fileID}_out-sorted.bam

#Calculate conversion rate, will comment out now to save time
#dnmtools bsrate [OPTIONS] -c <chroms> <input.sam>

echo '************************************'
echo Start ${fileID} methcounts:
echo '************************************'
date

mkdir $workdir/methcounts

#Calculating counts. Doing CpG Context only to save time
dnmtools counts -v -cpg-only -c $Index -o $workdir/methcounts/${fileID}.meth $workdir/uniqueBams/${fileID}_out-sorted.bam


echo '************************************'
echo Start ${fileID} Sym:
echo '************************************'
date

mkdir $workdir/symCpGs

#Output to symmetryic CpGs
dnmtools sym -o $workdir/symCpGs/${fileID}_sym.meth $workdir/methcounts/${fileID}.meth


echo '************************************'
echo Start ${fileID} Bigwigs:
echo '************************************'
date

mkdir $workdir/meth_bigwigs


#awk -v OFS="\t" '{print $1, $2, $2+1, $4":"$6, $5, $3}'  $workdir/symCpGs/${fileID}_sym.meth >  $workdir/symCpGs/${fileID}_symmeth.bed

cut -f 1-3,5 $workdir/symCpGs/${fileID}_symmeth.bed | wigToBigWig /dev/stdin /labs/vsebast/DJS/genomes/HomoSap/UCSC_ver/hg38.chrom.sizes $workdir/meth_bigwigs/${fileID}_meth.bw


echo '************************************'
echo Start ${fileID} LMRs:
echo '************************************'
date

mkdir $workdir/LMRs

#Finding HMRs/LMRs
dnmtools hmr -v -o $workdir/LMRs/${fileID}.hmr $workdir/symCpGs/${fileID}_sym.meth




######### making merged LMRs (separate to above script ran as an array in SCG) ###############

mkdir Merged_sym

######NEO Merged ########

dnmtools merge -o Merged_sym/NEO2_allPs.meth symCpGs/NEO2_P2_sym.meth symCpGs/NEO2_P5_sym.meth symCpGs/NEO2_P8_sym.meth

mkdir Merged_LMRs

#Finding HMRs/LMRs
dnmtools hmr -v -o Merged_LMRs/NEO2_allPs.hmr Merged_sym/NEO2_allPs.meth



#### old merged ###

dnmtools merge -o Merged_sym/OLD3_allPs.meth symCpGs/OLD3_P2_sym.meth symCpGs/OLD3_P5_sym.meth symCpGs/OLD3_P8_sym.meth

mkdir Merged_LMRs

#Finding HMRs/LMRs
dnmtools hmr -v -o Merged_LMRs/OLD3_allPs.hmr Merged_sym/OLD3_allPs.meth
