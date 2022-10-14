Mixedstrobes/Altstrobes
===========

This repository contains information, scripts, and sourcecode to run the analysis for various strobemers in Maier and Sahlin, 2022 [bioRxiv LINK TBD]. For info about strobemers, see [the repo](https://github.com/ksahlin/strobemers).


### What is a mixedstrobe?

Mixedstrobes is a mixed seeding technique that samples either a $k$-mer or a strobemer at a specified fraction. They may be sampled with all three strobemers seeding techniques (minstrobes, hybridstrobes and randstrobes) as well as altstrobes. Analogous to strobemers, we parameterize mixedstrobes as $(n,\ell,w_{min}, w_{max}, q)$, where $n$ is the number of strobes $n$, $\ell$ is the strobe length, $w_{min}$ and $w_{max}$ the minimum and maximum downstream offset to last window, and $q$ the strobemer fraction. Details on construction are found in [bioRxiv LINK TBD](XXX)

### What is an altstrobe?

Altstrobes are modified randstrobes where the strobe length is alternating between shorter and longer strobes. Contrary to randstrobes, where all strobes are of equal size (e.g. $k/2$ for strobemers of order 2), altstrobes are build up of one short strobe $k_s$ one longer strobe $k_l$, with $|k_s| + |k_l| = k$. We parameterize altstrobes as $(n,|k_s|/|k_l|, w_{min}, w_{max})$. Details on constructions are found in [bioRxiv LINK TBD](XXX)

### This repository

The repository consists of 

- functions to generate strobemers in C++
- functions to generate strobemers in Python
- a tool `StrobeMap` implemented in both C++ and Python
- a modified version of 'minimap2' for strobemers
- scripts used for the evaluations in the [paper LINK TBD](XXX)

# C++ functions

The C++ library `strobemers_cpp/index.[cpp/hpp]` contains functions for creating altstrobes, mixedstrobes, randstobes, hybridstrobes and minstrobes. 

You can copy the `index.cpp` and `index.hpp` files in the `strobemers_cpp` folder in this repository if you want to use either altstrobes, mixedstrobes, randstobes, hybridstrobes, or minstrobes in your project. They have the same user interface as described [here](https://github.com/ksahlin/strobemers#c-functions).


# Python functions

The `indexing_Maier_altstrobes.py` module located in the `modules` folder contains functions for generating k-mers, spaced k-mers, minimizers, and strobemers (minstrobes, hybridstrobes, randstrobes, mixedstrobes, and altstrobes) of any order. They have the same user interface as described [here](https://github.com/ksahlin/strobemers#python-functions).  

# StrobeMap (C++)


The tool `StrobeMap_Maier` is a program which roughly has the same interface as `MUMmer`. `StrobeMap_Maier` takes a reference and query file in multi-fasta or fastq format. It produces NAMs (Non-overlapping Approximate Matches) between the queries and references and outputs them in a format simular to nucmer/MUMmer. See [supplementary material](https://genome.cshlp.org/content/suppl/2021/10/19/gr.275648.121.DC1/Supplemental_Methods.pdf) Section A in the paper for definition of NAMs.

## Installation

`StrobeMap` (currently still without altstrobes and mixedstrobes) is available through [bioconda](https://bioconda.github.io/recipes/strobemap/README.html#package-strobemap). You can also acquire precompiled binaries for Linux from [here](https://github.com/benjamindominikmaier/mixedstrobes altstrobes/tree/main/src/strobemers_cpp/binaries). For example, for linux, simply do

```
wget /home/benji/Desktop/Altstrobe Mixedstrobe/src/cpp/StrobeMap/binaries/StrobeMap
chmod +x StrobeMap-0.0.2
./StrobeMap-0.0.2  # test program
```

If you want to compile from source, you need to have a newer `g++` and `zlib` installed. Then do the following:

```
git clone https://github.com/ksahlin/strobemers
cd strobemers/strobemers_cpp/

# Needs a newer g++ version. Tested with version 8 and upwards.

g++ -std=c++14 main.cpp index.cpp -lz -fopenmp -o StrobeMap -O3 -mavx2
```  

If zlib is not already installed on your system, it can be installed through, e.g., conda by

```
conda install -c anaconda zlib
```

If you dont have conda, download and install [here](https://zlib.net/). 


## Usage

```
$ ./StrobeMap 

StrobeMap VERSION 0.0.2

StrobeMap [options] <references.fasta> <queries.fast[a/q]>
options:
	-n INT number of strobes [2]
	-k INT strobe length, limited to 32 [20]
	-v INT strobe w_min offset [k+1]
	-w INT strobe w_max offset [70]
	-f INT strobe_fraction [60]
	-t INT number of threads [3]
	-o name of output tsv-file [output.tsv]
	-c Choice of protocol to use; kmers, minstrobes, hybridstrobes, altstrobes, randstrobes [randstrobes]. 
	-C UINT Mask (do not process) strobemer hits with count larger than C [1000]
	-L UINT Print at most L NAMs per query [1000]. Will print the NAMs with highest score S = n_strobemer_hits * query_span. 
	-S Sort output NAMs for each query based on score. Default is to sort first by ref ID, then by query coordinate, then by reference coordinate. 
	-s Split output into one file per thread and forward/reverse complement mappings. 
	   This option is used to generate format compatible with uLTRA long-read RNA aligner and requires 
	   option -o to be specified as a folder path to uLTRA output directory, e.g., -o /my/path/to/uLTRA_output/ 
```

```
# randstrobes (3,30,31,60)
StrobeMap -k 30 -n 3 -v 31 -w 60 -c randstrobes -o mapped.tsv  ref.fa query.fa
```


## Output

The output is a file `matches.tsv` in the output folder. You can se a custom outfile name with the parameter `--prefix`.
Output format is a tab separated file on the same format as MUMmer, with identical fields except the last one which is approximate reference sequence match length instead of what MUMmer produce:

```
>query_accession
ref_id  ref_pos query_pos   match_length_on_reference
```


# Minimap2-Strobemers

We implemented subsampled randstrobes, mixedstrobes, and altstrobe seeds in minimap2. A precompiled binary for Linux can be directly downloaded from [here](https://github.com/benjamindominikmaier/mixedstrobes altstrobes/tree/main/src/minimap-strobemers/minimap2).

If you want to compile from the source, you need to have a C compiler, GNU make and zlib development files installed. Then type make in the source code directory to compile. If you see compilation errors, try make sse2only=1 to disable SSE4 code, which will make minimap2-strobemer slightly slower.

All default minimap2 settings can be found in the original minimap implementation manual [here](https://github.com/lh3/minimap2). Our implementation adds the following settings:
```
 --mixedstrobes  sample mixedrandstrobes80 instead of k-mers
 --altstrobes    sample altstrobes (k,2k) instead of k-mers
 --randstrobes   sample randstrobes (k,k) instead of k-mers
    -i INT       minimal strobe offset [25].
    -j INT       maximial strobe offset [50].
```

# Evaluation Scripts

tbd


<!-- ## What is a NAM?

The aim is to output regions that _approximately match eachother_ on the reference and query (just like MEMs are exact matches between query and reference). The NAM regions can then be used to detect as candidate regions for alignment, clustering, or any other downstream analysis. 

The 'approximate' part is that there is a strobemer match, and the maximal is that we will merge matches that overlap on both query and reference sequence (if the order of the matches is the same on both the query and reference sequence). The 'maximal' _seems_ to be well defined if there are no repetitive matches in either the query or the reference. However, is is not trivial to compute if there are nested repeats in _both_ the query and reference. Currently, this is what is implemented:

For kmers, any two k-mer matches spanning positions `(q_1, q_1+k)` and `(q_2, q_2+k`) on the query and positions `(r_1, r_1+k)` and `(r_2, r_2+k)` on the reference where `q_1 <= q_2 <= q_1+k <= q_2+k` and `r_1 <= r_2 <= r_1+k <= r_2+k` are merged into one match of length `r_2+k - r_1`. Any chain of such overlapping matches are merged into one match. 


For strobemers, `StrobeMap` saves the positions for both the first and second strobe. Two strobemers with start positions `(q_1, q'_1)` and `(q_2, q'_2)` on the query and `(r_1, r'_1)` and `(r_2, r'_2)` on the reference with length `k` strobes _overlap_ if `q_1 <= q_2 <= q'_1 +k` and `r_1 <= r_2 <= r'_2+k`. If there is an overlap the two strobes are merged into one match of length `max(q'_1+k, q'_2 + k) - q_1`. Notice that because of the random length between the strobes, we can either have `q_1 <= q_2 <= q'_1 <= q'_2` or `q_1 <= q_2 <= q'_2 <= q'_1`, hence we need the `max` function. Any chain of such overlapping matches are merged into one match. 


The tool currently have a known bug of not being able to merge matches when there exist a repeat occuring at least twice within _both_ the query and reference sequence. In this case the matches may become fragmented, i.e., not merged into MAMs. -->
<!-- 
## Proof of concept

I aligned ONT cDNA reads (meadian error rate 7.0%) from [this synthetic RNA dataset](https://www.ebi.ac.uk/ena/browser/view/PRJEB34849) to SIRV transcripts [available here](https://github.com/ksahlin/strobemers/blob/main/data/sirv_transcripts.fasta) using minimap2 with parameters `-k 10 -w 1` providing very sensitive/accurate alignment. The ONT reads have been processed into full length reads using pychopper. Thus, ideally all reads should span and align to the full transcript. I selected 100 reads aligning to each SIRV transcript (primary alignment), and compared match coverage and number of hits between the 100 reads and their reference transcript using kmers **(k=30)** and strobemers **(n=2,k=15,w50)** giving the same subsequence length of 30nt each.

The aim is that the matching should provide candidate regions/sequences to perform exact alignment against. Match coverage and number of hits are two important features for sequence matching. Match coverage in this experiment is the fraction of reference sequence covered. Number of hits is the number of MAMs per read. We want the match coverage to be high in this experiment since we know that the reads align well to the their respective SIRV reference. However, together with a high coverage, we want the number of matches to be as low as possible (where 1 is best) in order for fast post prosessing/clustering of matches (a.k.a. seeds) and low disk space. 

Below I show the match coverage and number of hits for strobemers and kmers in this experiment separated for each of the 63 SIRVs with more than 100 primary alignment. The line shows the mean and the shaded area around the line is the standard deviation of the data (i.e., coverage/nr matches) for each SIRV.


![match coverage](data/plot_coverage.png)
![number of hits](data/plot_nr_hits.png)

The two above metrics could be studied from another angle, which is the match length normalized with the SIRV transcript length. The plot below shows the mean normalized match length for kmers and strobemers.

![norm match length](data/plot_normalized_match_length.png)

### Window placement matters

Above plots were produced with a second strobe produced from a window adjacent to the first kmer `k_1`, i.e., at offset positions `[0,50]` of the end of `k_1`. If we place the window in which we sample the second strobe a bit further donwstream (here I choose `[20,70]`), we get the following improved results where many of the strobemer MAMs cover the complete reference.

![match coverage](data/plot_coverage_w20_70.png)
![number of hits](data/plot_nr_hits_w20_70.png)
![norm match length](data/plot_normalized_match_length_w20_70.png) -->


CREDITS
----------------

Kristoffer Sahlin, Effective sequence similarity detection with strobemers, Genome Res. November 2021 31: 2080-2094; doi: https://doi.org/10.1101/gr.275648.121
Benjamin D. Maier and Kristoffer Sahlin, Entropy predicts fuzzy seed sensititivity, bioRxiv, October 2022; XXX
