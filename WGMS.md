
## Overview
Using DNMTools to analyze whole-genome methylation sequencing data

### [DNMTools](https://github.com/smithlabcode/dnmtools)

## General setup
Using [SCG](https://ondemand.scg.stanford.edu/) for storage and computing:

Using specific server
- From shell: ssh tmurty@smsh11dsu-srcf-d15-35.scg.stanford.edu
- Mahdi goes to same server each time (35)

Using sessions to have separate jobs running simultaneously
- tmux -- gives you a session and you can then switch between sessions
- tmux control b, release everything, and then press d -- exit session
- tmux ls -- see all the sessions; pwd and hostname to see that connected to server
- tmux a -t15 (this connects to session 15, for example)

## Installing DNMTools
1. DNMTools requires [HTSLib](https://www.biostars.org/p/328831/)
6. wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
7. tar -vxjf htslib-1.9.tar.bz2
8. cd htslib-1.9
2. route to home directory
3. wget https://github.com/smithlabcode/dnmtools/releases/download/v1.2.4/dnmtools-1.2.4.tar.gz
4. tar -zxvf dnmtools-1.2.4.tar.gz
5. cd dnmtools-1.2.4 && mkdir build && cd build
6. ./configure --prefix=/home/tmurty/htslib-1.9
   - what worked for Mahdi: ./configure --prefix=/home/moqri/htslib-1.9/test
8. make
9. install
10. cd back to DNAMtools
11. ISSUE: Where is "include" and "lib"?
    - ./configure CPPFLAGS='-I /home/moqri/dnmtools-1.2.4/htslib-1.9/test/include'
             LDFLAGS='-L /home/moqri/dnmtools-1.2.4/htslib-1.9/test/lib'

   - ../configure CPPFLAGS='-I /home/tmurty/htslib-1.9/headers' \
             LDFLAGS='-L /home/tmurty/htslib-1.9/lib'



