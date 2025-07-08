```
module add sratoolkit
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5945461/SRR5945461
sam-dump SRR5945461 > SRR5945461.sam #use bam if don't need to check the data
#in R
install.packages("devtools")
library(devtools)
install_github("broadinstitute/ichorCNA")
```
