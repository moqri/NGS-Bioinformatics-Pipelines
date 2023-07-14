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
