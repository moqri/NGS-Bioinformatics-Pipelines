
## Overview
Using DNMTools to analyze whole-genome methylation sequencing data

### [DNMTools GitHub](https://github.com/smithlabcode/dnmtools)
### [DNMTools Read the Docs](https://dnmtools.readthedocs.io/en/latest/)

Using [SCG](https://ondemand.scg.stanford.edu/) for storage and computing:
[Primer for SCG](https://github.com/nicolerg/resources/blob/master/scg_primer.md) with information about nodes, modules, etc.

## General setup with installing DNMTools and HTSLib

Using specific server
- From shell: ssh tmurty@smsh11dsu-srcf-d15-35.scg.stanford.edu
- Mahdi goes to same server each time (35)

Using sessions to have separate jobs running simultaneously

```tmux``` -- gives you a session and you can then switch between sessions

tmux control b, release everything, and then press d -- exit session

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


### Mapping
1. Mapping using abismal
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
