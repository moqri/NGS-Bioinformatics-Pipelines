## Overview
Using DNMTools to analyze whole-genome methylation sequencing data

### [DNMTools GitHub](https://github.com/smithlabcode/dnmtools)
### [DNMTools Read the Docs](https://dnmtools.readthedocs.io/en/latest/)

Using [SCG](https://ondemand.scg.stanford.edu/) for storage and computing:
[Primer for SCG](https://github.com/nicolerg/resources/blob/master/scg_primer.md) with information about nodes, modules, etc.

## General setup with installing DNMTools and HTSLib

Using specific server
- From shell: ```ssh tmurty@smsh11dsu-srcf-d15-35.scg.stanford.edu```
- Mahdi goes to same server each time (35)

Using sessions to have separate jobs running simultaneously

```tmux``` -- gives you a session and you can then switch between sessions

within a tmux session: control b, release everything, and then press d -- exit session

```tmux ls``` -- see all the sessions; pwd and hostname to see that connected to server

```tmux a -t15``` (this connects to session 15, for example)

## Setup using SCG modules
1. Navigate to directory within vsebast/shared: cd /oak/stanford/scg/lab_vsebast/shared/wgms
2. Want to load module for DNMTools using SCG
3. ```module spider dnmtools```
4. ```module add dnmtools```
   - advantage: the above method of installation is having issues, so this is a single step of loading a module
   - disadvantage: we cannot (re)install these pacakges as the most updated version (stuck with whatever SCG uses)
5. Test that dnmtools is properly working by typing ```dnmtools``` (no quotes)
   - would not work in a separate session
  
## Using DNMTools
### Create an index using abismal: [github](https://github.com/smithlabcode/abismal/releases) & [docs](https://dnmtools.readthedocs.io/en/latest/abismal/)
1. Indexing the genome
```
module load dnmtools
dnmtools abismalidx hg38.fa hg38.idx
```
  
### Steps before mapping
1. Decompress raw .fq.gz file
   
   ```gzip -d <file.gz>```
   - Faster to work with decompressed files, but they take more space
   - abismal can work with compressed files (in the process decompresses and then deletes)
3. ```head <file.fq>```
   - can see the content of the file
   - the first line is the location of the reads
   - the second line is the actual read: no C because there is conversion in methylation data (compared to seeing C in DNA data)
   - ```head <file.fq> -n2``` will give 2 rows
4. Removing adapters (trimming)
5. Quality control (if you are not sure of quality of sequencing itself)


### Mapping using abismal (Mahdi's approach)
   - General: ```$dnmtools abismal [OPTIONS] input.fq [input-r2.fq]```
   - ```$dnmtools abismal -i $ind "$f"_1_val_1.fq.gz "$f"_2_val_2.fq.gz -t $p -v```
            - ind is the genome index
            - "$f"_1_val_1.fq.gz "$f"_2_val_2.fq.gz are the two paired read files
   - ```f=d62_M38_CKDL230014191-1A_H5NYWDSX7_L2```
      - Defining f in such a way that it is the same for the paired reads 
   - ```module add trim-galore/0.6.7```
   - ```trim_galore --paired -q 0 --length 0 "$f"_1.fq.gz  "$f"_2.fq.gz```
      - does quality control, removes adapters (i.e., all recommended steps)
      - output: trimmed version of file + _report.txt
      - slow process, but can make it 20-40x faster by doing multi process using "j"
            - ```trim_galore --paired -q 0 --length 0 "$f"_1.fq.gz  "$f"_2.fq.gz -j 60```

         - ```cat file_trimming_report.txt```
             - To evaluate the trimming
             - should see that there are ~0% C because of methylation data
             - adapter information: Illumina TruSeq (confirm this is what we used)
             - ```cat d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_1.fq.gz_trimming_report.txt```
    
   - ```module add samtools```
      - If not done yet in this session:  ```module add dnmtools``` & ```cd /labs/vsebast/shared/wgms```
   - ```dnmtools abismal -i <index> <val1> <val2> -t <nodes> -v | samtools view -b > mapped.bam```
       - ```dnmtools abismal -i hg38.idx d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_1_val_1.fq.gz d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_2_val_2.fq.gz -t 60 -v | samtools view -b > d62_M38_val_mapped.bam```
       - since the file is >10GB, then would expect the mapping to take ~hours (13GB and 14GB)
       - started at 3:44PM 07/14/2023; 2% completed at 3:50PM; 10% completed at 4:18pm; 20% at 4:54 -- approximate 6 hours total
         ```
         (base) [tmurty@smsh11dsu-srcf-d15-35 wgms]$ dnmtools abismal -i hg38.idx d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_1_val_1.fq.gz d62_M38_CKDL230014191-
         1A_H5NYWDSX7_L2_2_val_2.fq.gz -t 60 -v | samtools view -b > d62_val_mapped.bam                     
         [Fri Jul 14 15:43:41 2023] using 60 threads to map reads.
         [Fri Jul 14 15:43:41 2023] input (PE): d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_1_val_1.fq.gz, d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_2_val_2.fq.gz
         [Fri Jul 14 15:43:41 2023] output (SAM): [stdout]
         [Fri Jul 14 15:43:41 2023] loading index hg38.idx
         [Fri Jul 14 15:43:50 2023] loading time: 8.818743s
         [mapping reads|===================================================|100%]
         [Fri Jul 14 21:25:36 2023] total mapping time: 20506.728474s
         ```
### Mapping using abismal (Readthedocs) - July 15 TM approach
1. For paired-end reads:
   - ```abismal -i <index> -o <output SAM> [options] <input fastq 1> <input fastq 2>```
   - ``` dnmtools abismal -i hg38.idx -o d62_M38_mapped.sam -t 60 -v d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_1_val_1.fq.gz d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_2_val_2.fq.gz ```
   - <4 hours to map to .sam
     ```
     (base) [tmurty@smsh11dsu-srcf-d15-35 wgms]$ dnmtools abismal -i hg38.idx -o d62_M38_mapped.sam -t 60 -v d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_1_val_1.fq.gz       d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_2_val_2.fq.gz
      [Sat Jul 15 16:15:34 2023] using 60 threads to map reads.
      [Sat Jul 15 16:15:34 2023] input (PE): d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_1_val_1.fq.gz, d62_M38_CKDL230014191-1A_H5NYWDSX7_L2_2_val_2.fq.gz
      [Sat Jul 15 16:15:34 2023] output (SAM): d62_M38_mapped.sam
      [Sat Jul 15 16:15:34 2023] loading index hg38.idx
      [Sat Jul 15 16:15:42 2023] loading time: 8.277908s
      [mapping reads|===================================================|100%]
      [Sat Jul 15 19:36:51 2023] total mapping time: 12068.419744s 
      ```

2. Formatting [ReadtheDocs](https://dnmtools.readthedocs.io/en/latest/format/)
   - ```dnmtools format -f abismal input.bam output.sam```
      - For Mahdi's approach mapped bam: d62_M38_val_mapped.bam
      - ```dnmtools format -f abismal d62_M38_val_mapped.bam d62_M38_val_mapped_format.sam ```
   - ```dnmtools format -f abismal -t 8 -B input.bam output.bam```

### Calculating DNAm levels
   - ```counts``` which was previuosly ```methcounts```
   - ```$ dnmtools counts -c /path/to/genome.fa -o output.meth input.sam```
      - can input the .sam file (can be read as text file, relatively large) or compressed .bam file (binary)
      - "The input mapped reads file (input.bam) is in SAM/BAM format. The reads should be sorted so those mapping to the same chromosome are consecutive in the file. Duplicate reads should be probably be removed first, but that depends on your data." [DNMTools](https://dnmtools.readthedocs.io/en/latest/counts/)
      - ```$ dnmtools counts -c /labs/vsebast/shared/wgms/hg38.fa -o d62_M38.meth d62_M38_val_mapped.bam``` OR ``` dnmtools counts -c /labs/vsebast/shared/wgms/hg38.fa -o d62_M38_fromsam.meth d62_M38_mapped.sam ``` 
         - Get out: ```reads in SAM file not sorted ```
         - **Tara's approach to sort**: attempted solution to be able to run counts based on error above [Biostars forum](https://www.biostars.org/p/319730/)
            - ``` samtools sort d62_M38_val_mapped.bam -m 1G -@ 8 -o d62_M38_val_mapped_sorted.bam```
               - Need to specify -m 1G otherwise each temp sort file is 145 MB and there are many many created
               - will see this: ```[bam_sort_core] merging from 15 files and 8 in-memory blocks...```
               - output file is 16GB vs pre-sorted (mapped) is 25GB
            - Using sorted file for counts
               - ``` dnmtools counts -c /labs/vsebast/shared/wgms/hg38.fa -o d62_M38.meth d62_M38_val_mapped_sorted.bam ```
               - results in 36G file ```d62_M38.meth ```

### Global Methylation Summary Stats (Tara following ReadtheDocs 7/17/23)
- [Levels](https://dnmtools.readthedocs.io/en/latest/levels/)
- ```dnmtools levels -o output.levels input.meth```
- ```dnmtools levels -o d62_M38.levels d62_M38.meth```
- ``` view d62_M38.levels ```

### Hypomethylated Regions (HMR) (Tara following ReadtheDocs 7/17/23)
- [hmr](https://dnmtools.readthedocs.io/en/latest/hmr/)
- "Running hmr requires a file of methylation levels formatted like the output of the counts. For identifying HMRs in mammalian methylomes, use the symmetric CpG methylation levels. This is obtained by using the sym command after having used the counts command."
   - [sym info](https://dnmtools.readthedocs.io/en/latest/sym/)
      - "The above command will merge all CpG pairs while also discarding sites with an indication that the CpG has mutated. Note that as long as one site of the pair is mutated, the pair is discarded. This is the default mode." How should mutations be taken into account? 
   - ```dnmtools sym -o human_esc_CpG.meth human_esc.meth```
   - ```dnmtools sym -o d62_M38_CpG.meth d62_M38.meth``` where _CpG.meth is version with collapsed counts for symmetric CpGs sites
      - results in 863M output file
- ```dnmtools hmr -p params.txt -o output.hmr input.meth```
- ```dnmtools hmr -p params_d62_M38.txt -o d62_M38.hmr d62_M38_CpG.meth```
   - using d62_M38.meth (output from counts) gives error: ```error: input is not symmetric-CpGs: d62_M38.meth```
   - d62_M38.hmr is 1.4M 
