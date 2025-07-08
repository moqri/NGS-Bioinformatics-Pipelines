```
module add sratoolkit
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5945461/SRR5945461
sam-dump SRR5945461 > SRR5945461.sam
#in R
install.packages("devtools")
library(devtools)
install_github("broadinstitute/ichorCNA")
```
