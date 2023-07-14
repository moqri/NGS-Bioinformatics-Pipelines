## ~~Installing DNMTools~~ (not completed because of issues -- used SCG module)
1. DNMTools requires [HTSLib](https://www.biostars.org/p/328831/)
6. wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
7. tar -vxjf htslib-1.9.tar.bz2
8. cd htslib-1.9
2. route to home directory
3. wget https://github.com/smithlabcode/dnmtools/releases/download/v1.2.4/dnmtools-1.2.4.tar.gz
4. tar -zxvf dnmtools-1.2.4.tar.gz
5. cd dnmtools-1.2.4 && mkdir build && cd build
6. cd /home/tmurty/htslib-1.9
7. ./configure --prefix=/home/tmurty/htslib-1.9
   - What worked for Mahdi, but did not work for Tara: ./configure --prefix=/home/moqri/htslib-1.9/test
8. make
9. install
10. cd back to DNAMtools
11. **ISSUE**: Where is "include" and "lib"?
    - ./configure CPPFLAGS='-I /home/moqri/dnmtools-1.2.4/htslib-1.9/test/include'
             LDFLAGS='-L /home/moqri/dnmtools-1.2.4/htslib-1.9/test/lib'

   - ../configure CPPFLAGS='-I /home/tmurty/htslib-1.9/headers' \
             LDFLAGS='-L /home/tmurty/htslib-1.9/lib'

### Calculating DNAm levels
   - ```counts``` which was previuosly ```methcounts```
   - ```$ dnmtools counts -c /path/to/genome.fa -o output.meth input.sam```
      - can input the .sam file (can be read as text file, relatively large) or compressed .bam file (binary)
      - "The input mapped reads file (input.bam) is in SAM/BAM format. The reads should be sorted so those mapping to the same chromosome are consecutive in the file. Duplicate reads should be probably be removed first, but that depends on your data." [DNMTools](https://dnmtools.readthedocs.io/en/latest/counts/)
   - using a previously generated bam file to test (since mapping is ongoing from previous section)
      - ```/labs/vsebast/moqri/data/PRJEB28044/60_Hmp01_blood_young.bam```
   - ```dnmtools counts -c hg38.fa -o DNAm.meth /labs/vsebast/moqri/data/PRJEB28044/60_Hmp01_blood_young.bam```
      - ~~unfortunately this does not run because it does not find chromosome 1 (likely it should be chr1 and it is just 1)~~
