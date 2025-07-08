```
module add sratoolkit
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5945461/SRR5945461
sam-dump SRR5945461 > SRR5945461.sam #use bam if don't need to check the data
#in R
install.packages("devtools")
library(devtools)
install_github("broadinstitute/ichorCNA")

 bitops, RCurl, HMMcopy

 install.packages("HMMcopy", lib = "lib")


sort -k1,1 -k2,2n SRR5945461.bedgraph > SRR5945461.sorted.bedgraph

awk  'BEGIN{chr=""} {if(chr!=$1){chr=$1; print "variableStep chrom="chr} print $2+1"\t"$4}' SRR5945461.sorted.bedgraph > SRR5945461.wig


```
