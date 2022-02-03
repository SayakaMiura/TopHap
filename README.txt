# TopHap_v1.2.1  

(Copyright 2022, Authors and Temple University; see license below)
Updated February 3, 2022
==================

## Description
TopHap infers bootstrap-supported phylogenies of common haplotypes in the given data. See Caraballo et al. (ref. 1) for the detail. The TopHap program has been developed by Sudhir Kumar. It is written in Python. You are free to download, modify, and expand this code under a permissive license similar to the BSD 2-Clause License (see below).

## Dependencies
1. python 3 (v3.7.9 and v3.8.3 were tested)
 python packages: 
    numpy
    biopython
 Note: If the installation of these python packages is not easy, you may want to use Anaconda for Python 3 (https://www.anaconda.com/distribution/). Or you can try python3 -m pip install [package name].

2. R (v3.5.3 and v4.0.3 were tested)
 R package: 
    ape
    phangorn 
Please make sure “Rscript” command is functional.

3. MEGA
 Please download the latest version from https://www.megasoftware.net/.

4. RaxML (optional)
 It can be downloaded at https://cme.h-its.org/exelixis/web/software/raxml/.

##How to use TopHap
1. Run vcf_json_parse.py. (optional)

Given the alignment of all genomes (i.e., aligned with Wuhan1 reference genome sequence), common nucleotides (positions with desired minor allele frequency (maf) threshold (e.g., > 5%) are extracted and haplotype alignments that contain only these genomic positions are generated for each spatiotemporal slice of the dataset (country and sampling month).  
python3 vcf_json_parse.py [input full genome alignment] --reference Wuhan1Gnome.fasta --min_subgroup_size 500 --one_based --skip_mismatches --min_freq 0.05 -o [output directory]

 ====input full genome alignment====
 Each genome sequence information needs to be saved in a json format. The format should be,
 {"variants":[[genomic position,nucleotide,count,frequency],...],
 "sequences":[{"V":[genomic position,nucleotide,...],"I":[sampling time,continent,country,PANGO ID,NextStrain ID]},...]},
 }
 For example,
 {"variants":[[13467,"-",1106856,0.9999936758207655],[23403,"G",1091067,0.9857290378303367]...],
 "sequences":[{"V":[1059,"T",3037,"T"],"I":["20200308","Africa","Algeria","B.1","20C"]},{"V":[10582,"T",13467,"-",14408,"T",23403,"G",25563,"T"],"I":["20200302","Africa","Algeria","B.1","20C"]},...]
 }
 Please do not use space and symbols for the country name. An example JSON file is provided, "Example.json." 

 ====options====
--reference [reference genome sequence file]: Reference genome sequence to determine mutant and non-mutant nucleotides.   
--min_subgroup_size [number]: The minimum genome count per slice (alignment file). When the count in a monthly slice is smaller than this minimum count, the slice is merged to an adjacent month. When the count is still smaller than the minimum count after merging slices, such slices will be ignored to compute minor allele frequency (maf).  
--one_based: If the genomic position in the JSON file is counted from 1, use this option.
--skip_mismatches: If gaps (-) are desired to be ignored for the calculation of maf, use this option.
--min_freq [a value between 0 and 1]: Minor allele frequency (maf) threshold. 
-o [output directory]: The directory all output files are stored.

 ====output files====
 (a) fasta files
 An alignment of haplotypes is produced for each spatiotemporal slice.
 (b) Haplotypes.txt
 Genomic positions are used to generate each fasta file. The genomic positions are counted from 0, and the order of the positions is the same as those in the fasta file.

 ====example====
 "Example.json" can be processed by,
python3 vcf_json_parse.py Example.json --reference Wuhan1Gnome.fasta --min_subgroup_size 500 --one_based --skip_mismatches --min_freq 0.05 -o ExampleHap 
 The output alignment files together with genomic positions extracted ("Haplotypes.txt") are stored in "ExampleHap" directory.

2. Run TopHap.py.
The main program that infers bootstrap-supported phylogenies of common haplotypes in the given data.
python3 TopHap.py [haplotype frequency cutoff] [number of bootstrap replicates] –Hap [path to the directory of the haplotype alignments] 

 ====options====
haplotype frequency cutoff: haplotypes with the desired frequency thresholds are selected.   
the number of bootstrap replicates: The desired number of bootstrap samples are generated.

 ====output files====
 (a) Bootstrap scored phylogenies (TopHap_bootstrap*.nwk)
 If there is more than one equally parsimonious tree, each tree is scored. Please select the one with the highest support.
 (b) Haplotype alignment (TopHap_prune.fasta)
 It can be used together with a Bootstrap scored phylogeny to infer ancestral states, from which you can find the timing of each mutation. The ancestral reconstruction function in MEGA (https://www.megasoftware.net/) is useful for this purpose. 
 (c) Results of each bootstrap replicate dataset are found in the Bootstrap directory.

 ====example====
 Example datasets can be found in Alignment.
 	python3 TopHap.py 0.05 100 -Hap Alignment
 
3. Run TopHap_Attach.py (optional)
Attach minor haplotype sequences into a TopHap phylogeny. For a TopHap phylogeny, haplotype of interest will be attachedpython3 TopHap_Attach.py [TopHap alignment] [TopHap tree] [minor haplotype] [path to raxml]

 ====inputs====
 TopHap alignment: A fasta file that contains TopHap haplotypes
 TopHap tree: A nwk file that contains TopHap tree
 minor haplotype: A fasta file that contains minor haplotypes that are attached to the TopHap tree

 ====output files====
 (a) A phylogeny with a minor haplotype attached (RaxMLoutputs\*.nwk)
 All trees are found at RaxMLoutputs directory.
 (b) RaxML output files
 Direct output files from RaxML are stored in RaxMLoutputsAll. 
 
 ====example====
 Example datasets can be found in Example_attach. TopHap alignment and tree are TopHap_prune.fasta and TopHap_bootstrap1.nwk, respectively. To attach haplotypes in MinorHap.fasta, run,
python3 TopHap_Attach.py Example_attach\TopHap_prune.fasta Example_attach\TopHap_bootstrap1.nwk Example_attach\MinorHap.fasta raxmlHPC-AVX.exe

===================
Reference:
[1] Marcos A. Caraballo-Ortiz, Sayaka Miura, Sergei L. K. Pond, Qiqing Tao, and Sudhir Kumar. TopHap: TopHap: Rapid inference of key phylogenetic structures from common haplotypes in large genome collections with limited diversity (2021) Submitted to Bioinformatics

--------
Copyright 2022, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
