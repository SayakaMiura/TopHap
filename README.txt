TopHap_v0.1.1  
(Copyright 2021, Authors and Temple University; see license below)

Updated January 24, 2021
==================

The TopHap program has been developed by Sayaka Miura. It is written in Python3.7 for Windows. You are free to download, modify, and expand this code under a permissive license similar to the BSD 2-Clause License (see below). 
TopHap identifies common haplotypes from alignments of SARS-CoV-2 genome sequences and infer bootstrap-supported phylogenies. Users can expand the phylogeny by providing additional alignments. See Caraballo et al. (ref. 1) for the detail. Sequences sampled from late December 2019 to early October 2020 are already included. So, TopHap will identify common haplotypes in alignments that are provided by users and infer bootstrap-supported phylogenies with the additional haplotypes.  

Dependencies
==================
1. python 3 (v3.7.9 and v3.8.3 were tested)
 python packages: 
    numpy
    biopython
 Note: If the installation of these python packages is not easy, you may want to use Anaconda for python 3.7 (https://www.anaconda.com/distribution/). Or, you can try python3-pip.

2. R (v3.5.3 and v4.0.3 were tested)
 R package: 
    ape
 Please make sure “Rscript ” command is functional.

3. MEGA
 Please download the latest version from https://www.megasoftware.net/.

How to use 
==================
1. Unzip Alignment.zip
 Make sure the TopHap main directory has Alignment directory that contains fasta files and text files.

2. Align your SARS-CoV-2 sequences with the NCBI reference sequence
 Place your alignment files in the directory of Alignment. 

3. Add the information of your alignment files in AdditionalData.txt 
 As an example, we added two additional alignments, named SA501.fasta and UK501.fasta. Your fasta file name should be listed in the first column and you can add a note in the second column.

4. (optional) List genomic positions that you like to include to construct haplotypes in MinorVariants.txt.
 As an example, we listed positions that were reported among UK sequences with 501 variant.

5. (optional) List genomic positions that you like to exclude to construct haplotypes in Excluded.txt.
 We already listed genomic positions that are excluded in our database, i.e., non-coding regions and ORF7. You can list additional positions in the same file. 

6. (optional) Change parameter setting.
 The default setting is listed below. If you wish to use different setting please edit lines 5-8 in TopHap.py, accordingly.
 vf (variant frequency cutoff): 0.05 (position with >5% frequencies are used)
 hf (haplotype frequency cutoff): 0.05 (haplotypes with >5% frequencies are selected)   
 bootstrap replicates: 100 (100 bootstrap samples are generated)
 hf for bootstrap replicates: 0.05 (for bootstrap sample, haplotypes with >5% frequencies are selected) 

7. Run TopHap.py.
 Python3 TopHap.py 

Output file
==================
The output files are found in the same directory as TopHap. 
1. Bootstrap scored phylogenies (TopHap_bootstrap*.nwk)
 If there are more than one equally parsimonious tree, each tree is scored. Please select one with the highest supports.

2. Haplotype alignment (TopHap_prune.fasta)
 It can be used together with a Bootstrap scored phylogeny to infer ancestral states, from which you can find the timing of each mutation. The ancestral reconstruction function in MEGA (https://www.megasoftware.net/) is useful for this purpose. 

3. Summary tables
 3.1 CommonVarPosLs.txt
 List genomic positions with >5% variant frequencies in subdata. The order of the positions in this table is same as the order in the Haplotype alignment (TopHap_prune.fasta). Also, subdata with >5% frequencies were listed.
 Note for subdata code (sequences with the same sampling month):
 T1: Dec 2019 and Jan 2020
 T2: Feb 2020
 T3: March 2020
 T4: April 2020
 T5: May 2020
 T6: June 2020
 T7: July 2020
 T8: Aug 2020
 T9: Sep and early Oct 2020

 3.2 TopHap.txt
 For each common haplotype, subdata IDs with >5% frequencies are listed. 

 3.3 TopHap_anno.txt
 For each input sequence, haplotype ID is given. When none of common haplotype was matched, "NA" is given.

 3.4 TopHap_count.txt
 For each common haplotype, the count of sequences is given.

==================
Reference:
[1] Marcos A. Caraballo-Ortiz, Sayaka Miura, Sergei L. K. Pond, Qiqing Tao, and Sudhir Kumar. TopHap: Building strain phylogenies using major haplotypes in large genome collections (2021) BioRxiv

--------
Copyright 2021, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

