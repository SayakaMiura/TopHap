TopHap_v1.1.1  
(Copyright 2021, Authors and Temple University; see license below)

Updated November 16, 2021
==================

TopHap infers bootstrap-supported phylogenies of common haplotypes in the given data. See Caraballo et al. (ref. 1) for the detail. The TopHap program has been developed by Sayaka Miura. It is written in Python3.7 for Windows. You are free to download, modify, and expand this code under a permissive license similar to the BSD 2-Clause License (see below). 

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
 Please make sure “Rscript” command is functional.

3. MEGA
 Please download the latest version from https://www.megasoftware.net/.

How to use 
==================
1. Align all genomes with outgroup sequences.

2. Convert the full alignment into an alignment of haplotypes.
 When information on sampling location and time of haplotypes is available, please select common variants (positions with desired minor allele frequency threshold, e.g., > 5%) for each spatiotemporal slice of the dataset that is regionally (e.g., continent, country, or city) and/or temporally (e.g., monthly) partitioned. Then, pool all the variant positions and generate an alignment of haplotypes for each spatiotemporal slice. All alignments need to include the same genomic positions with the same order. For each alignment, please use the FASTA format and name it as “*_Hap.fasta,” e.g., slice1_Hap.fasta. Create a new directory and place all the files in that directory. Please do not place unnecessary other FASTA files in this directory. Example datasets can be found in Alignment.
 
3. Prepare an alignment of haplotypes of outgroup sequences. 
 Please use the same genomic positions with the same order as the other haplotype alignments prepared above. Please name the file as “OutG.fasta” and place it in the same directory as the other alignments. Example OutG.fasta can be found in Alignment.


6. (optional) Change parameter setting.
 The default setting is listed below. If you wish to use different setting, please edit lines 8-10 in TopHap.py, accordingly.
 hf (haplotype frequency cutoff): 0.05 (haplotypes with >5% frequencies are selected)   
 bootstrap replicates: 100 (100 bootstrap samples are generated)
 hf for bootstrap replicates: 0.05 (for bootstrap sample, haplotypes with >5% frequencies are selected) 

7. Run TopHap.py.
 python3 TopHap.py –Hap [path to the directory of the haplotype alignments]
 Example datasets can be found in Alignment. To run,
 python3 TopHap.py -Hap Alignment

Output file
==================
The output files are found in the same directory as TopHap. 
1. Bootstrap scored phylogenies (TopHap_bootstrap*.nwk)
 If there are more than one equally parsimonious tree, each tree is scored. Please select one with the highest supports.

2. Haplotype alignment (TopHap_prune.fasta)
 It can be used together with a Bootstrap scored phylogeny to infer ancestral states, from which you can find the timing of each mutation. The ancestral reconstruction function in MEGA (https://www.megasoftware.net/) is useful for this purpose. 

3. Results of each bootstrap replicate dataset are found in Bootstrap directory.

===================
Reference:
[1] Marcos A. Caraballo-Ortiz, Sayaka Miura, Sergei L. K. Pond, Qiqing Tao, and Sudhir Kumar. TopHap: TopHap: Rapid inference of key phylogenetic structures from common haplotypes in large genome collections with limited diversity (2021) Submitted to Bioinformatics

--------
Copyright 2021, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
