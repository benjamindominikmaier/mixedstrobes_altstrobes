Mixedstrobes/Altstrobes
===========

### What is a strobemer? 

A strobemer is a seed extracted from text and used for text similarity searches described [here](https://genome.cshlp.org/content/31/11/2080). Strobemers are _fuzzy_ seeds and can match over mutations in text. They were initially designed for biological sequence analysis where point mutations and small shifts (insertions/deletions) are common, but they are agnostic to the alphabet of text. 

A small illustration below of two biological sequences `T` and `T'` with substitutions and insertions/deletions between them, and four strobemer seeds (`s_1`-`s_4` and `s'_1`-`s'_4`) extracted from each sequence. Mutations destroy some seeds (`s_1`/`s'_1`), but some remain intact, allowing retrieval of similar regions. The seeds have total length of 24 letters (two 12-letter strings separated by a pseudorandom spacing), but no consecutive 24 letter string (k-mer) is found between `T` and `T'`.

![strobe](https://user-images.githubusercontent.com/1714667/149466166-1f558539-730a-4b1b-9876-a28f24bb01a2.png)

### What is a mixedstrobe?

Mixedstrobes is a mixed seeding technique that samples either a $k$-mer or a strobemer at a specified fraction. They may be sampled with all three strobemers seeding techniques (minstrobes, hybridstrobes and randstrobes) as well as altstrobes. Analogous to strobemers, we parameterize mixedstrobes as $(n,\ell,w_{min}, w_{max}, q)$, where $n$ is the number of strobes $n$, $\ell$ is the strobe length, $w_{min}$ and $w_{max}$ the minimum and maximum downstream offset to last window, and $q$ the strobemer fraction.

The decision about whether a strobemer or a $k$-mer is seeded depends on the hash value of the first strobe $h(S[i:i+\ell])$ and a user-defined strobe fraction $q$. A pseudocode to construct mixedstrobes is given in the Suppl. section S1 of Maier & Sahlin (2022).

### What is an altstrobe?

Altstrobes are modified randstrobes where the strobe length is alternating between shorter and longer strobes. Contrary to randstrobes, where all strobes are of equal size (e.g. $k/2$ for strobemers of order 2), altstrobes are build up of one short strobe $k_s$ one longer strobe $k_l$, with $|k_s| + |k_l| = k$. We parameterize altstrobes as $(n,|k_s|,|k_l|, w_{min}, w_{max})$. We refer to sampled altstrobe with $n=2$ as $(|k_s|,|k_l|)$ or $(|k_l|,|k_s|)$, depending on if the short strobe was used first or second, respectively. Whether the first strobe is of length $|k_s|$ or $|k_l|$ is decided based on the hash value of the substring of length $|k_s|$ (\textit{i.e.}, the potential first strobe). Altstrobes are implemented for all even orders. In the default version, long strobes $kl$ have length $2*k_s$. A pseudocode to construct altstrobes is given in the Suppl. section S1 of Maier & Sahlin (2022).

### This repository

The repository consists of 

- functions to generate strobemers in C++
- functions to generate strobemers in Python
- a tool `StrobeMap` implemented in both C++ and Python
- a modified version of 'minimap2' for strobemers
- scripts used for the evaluations in the [paper](XXX)


### Other implementations of strobemers 

Other strobemer implementations are found here
- [C++](https://github.com/BGI-Qingdao/strobemer_cpptest)
- [Go](https://github.com/shenwei356/strobemers)

The construction time is dependent on the implementation. The times reported in the preprint are for the C++/Python implementations in this repository.

# C++ functions

The C++ library `strobemers_cpp/index.[cpp/hpp]` contains functions for creating altstrobes, mixedstrobes, randstobes, hybridstrobes and minstrobes. 

You can copy the `index.cpp` and `index.hpp` files in the `strobemers_cpp` folder in this repository if you want to use either altstrobes, mixedstrobes, randstobes, hybridstrobes, or minstrobes in your project. The implementation of these functions uses bitpacking and some other clever tricks (inspired by [this repo](https://github.com/lh3/kmer-cnt)) to be fast. Because of bitpacking, the limitation is that any single strobe cannot be lager than 32, which means that the **maximum strobemer length allowed in this implementation** is `3*32 = 96` for strobemers of order 3, and `2*32 = 64` for order 2. This should be large enough for most applications. 

The functions in the `index.cpp` file can be used as follows:

```
#include "index.hpp"

typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> strobes_vector;
strobes_vector randstrobes3; // (kmer hash value, seq_id, strobe1_pos, strobe2_pos, strobe3_pos)
seq = "ACGCGTACGAATCACGCCGGGTGTGTGTGATCGGGGCTATCAGCTACGTACTATGCTAGCTACGGACGGCGATTTTTTTTCATATCGTACGCTAGCTAGCTAGCTGCGATCGATTCG";
n=3;
k=15;
w_min=16;
w_max=30;
seq_id = 0; // using integers for compactness, you can store a vector with accessions v = [acc_chr1, acc_chr2,...] then seq_id = 0 means v[0].

randstrobes3 = seq_to_randstrobes3(n, k, w_min, w_max, seq, seq_id);
for (auto &t : randstrobes3) // iterate over the strobemer tuples
{
strobemer_hash = std::get<0>(t);
strobe1_pos = std::get<2>(t);
strobe2_pos = std::get<3>(t);
strobe3_pos = std::get<4>(t);
// if you want the actual strobemer sequences:
randstrobe = seq.substr(strobe1_pos, k) + seq.substr(strobe2_pos, k)+ seq.substr(strobe3_pos, k);
}
```

If you are using some of `seq_to_randstrobes2`, `seq_to_hybridstrobes2`, or `seq_to_minstrobes3` they return the same vector tuples but position of strobe 2 copied twice, i.e., `(kmer hash value, seq_id, strobe1_pos, strobe2_pos, strobe2_pos)`. There is no reason for this and for any high performance application the function could be modified to return the minimal needed information.

The benchmarking in Table S3 in the [supplemental methods](https://genome.cshlp.org/content/suppl/2021/10/19/gr.275648.121.DC1/Supplemental_Methods.pdf) found that randstrobes is roughly as fast as hybridstrobes and minstrobes. Furthermore, randstrobes is unexpectedly fast in this implementation in general, about 1.5-3 times slower than generating k-mers for randstrobes of (n=2, s=15, w_min=16,w_max=70). What takes time is pushing the tuples to a vector, not computing the strobemers. Bare construction time could be further compared if preallocating an array of fixed size to remove the resizing when pushing to vectors. Nevertheless, the takehome message is that the generation of strobemers could be implemented so that it is fast, and will likely not be a bottleneck in most algorithms using them.

#### Notes for sequence mapping

The preprint describes shrinking the sampling windows `[w_min, w_max]` at ends of sequences to assure that a similar number of strobemers and k-mers created. However, in, e.g., read mapping, there is little to no gain in shrinking windows. This is because if we shrink windows at the ends of reads, the strobemer extracted from the read in those windows cannot be guaranteed to (but may) be the same as in the reference, as first described in [this issue](https://github.com/ksahlin/strobemers/issues/2). The more the window(s) are shrunk, the less likely the strobers are to match between the sequences, and the probability of matching a strobemer after the last window (original size) is completely outside the read is 0 (if disregarding false matches). After noting this, my implementation only shrink the last strobemer window regardless of the number of strobes (i.e., there is a positive probablility of a match even if window is shrunk). This means that there will be `n - (k + w_min) +1` strobemers of order 2 generated form a sequence of length `n`, and `n - (k + w_max + w_min) +1` strobemers of order 3. In other words, we will only slide last strobe's window outside the sequence. Once it is fully outside the sequence we stop (illustrated in approach B for order 2 in [here](https://github.com/ksahlin/strobemers/issues/2).


# Python functions

The `indexing_Maier_altstrobes.py` module located in the `modules` folder contains functions for generating k-mers, spaced k-mers, minimizers, and strobemers (minstrobes, hybridstrobes, and randstrobes) of any order. For randstrobes, there are two ways to create them. The first way is with the function `randstrobes`, which takes a string, k-mer size, and upper and lower window limits and returns a dictionary with positions of the strobes as keys and the hash value of the randstrobe sequence (strings) as values. For example

```
from modules import indexing
all_mers = defaultdict(list)
for (p1,p2,p3), h in indexing.randstrobes(seq, k_size, w_min, w_max, order = 3).items():
    # all_mers is a dictionary with hash values as keys and 
    # a list with position-tuples of where the strobemer is sampled from
    all_mers[h].append( (p1,p2,p3) )  
```
Functions `minstrobes`, `hybridstrobes` and `altstrobes`  have the same interface.

The second way is to call `randstrobes_iter` which is a generator. Similarly to `randstrobes`, `randstrobes_iter` takes a string, k-mer size, and upper and lower window size, but instead yields randstrobes from the sequence and is not as memmory requiring as the `randstrobes` function which store and returns all the strobes in a dictionary. `randstrobes_iter` generating randmers of order 2 can be used as follows

```
from modules import indexing
all_mers = defaultdict(list)
for (p1,p2), s in indexing.randstrobes_iter(seq, k_size, w_min, w_max, order = 2, buffer_size = 1000000):
    all_mers[s].append( (p1,p2) )  
```
Functions `minstrobes_iter`, `hybridstrobes_iter`, `altstrobes_iter` as well as their corresponding mixedstrobes have the same interface.

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

### Common installation from source errors

If you have `zlib` installed, and the `zlib.h` file is in folder `/path/to/zlib/include` and the `libz.so` file in `/path/to/zlib/lib` but you get 

```
main.cpp:12:10: fatal error: zlib.h: No such file or directory
 #include <zlib.h>
          ^~~~~~~~
compilation terminated.
```

add `-I/path/to/zlib/include -L/path/to/zlib/lib` to the compilation, that is

```
g++ -std=c++14 -I/path/to/zlib/include -L/path/to/zlib/lib main.cpp index.cpp -lz -fopenmp -o StrobeMap -O3 -mavx2
``` 


# StrobeMap (Python)

`StrobeMap` implements order 2 and 3 randstrobes (default), minstrobes, hybridstrobes, altstrobes as well as kmers. The tool produces NAMs (Non-overlapping Approximate Matches; see explanation in preprint) for both strobemers and kmers.
Here are some example uses:

```
# Generate hybridstrobe matches (hybridstrobe parametrization (2,15,20,70)) 
# between ONT SIRV reads and the true reference sequences

./StrobeMap --queries data/sirv_transcripts.fasta \
           --references data/ONT_sirv_cDNA_seqs.fasta \
           --outfolder strobemer_output/  --k 15 
           --strobe_w_min_offset 20 --strobe_w_max_offset 70


# Generate kmer matches (k=30) 
# between ONT SIRV reads and the true reference sequences

./StrobeMap --queries data/sirv_transcripts.fasta \
           --references data/ONT_sirv_cDNA_seqs.fasta \
           --outfolder kmer_output/  --k 30 --kmer_index

# Reads vs reads matching using randstrobes

./StrobeMap --queries data/ONT_sirv_cDNA_seqs.fasta \
           --references data/ONT_sirv_cDNA_seqs.fasta \
           --outfolder strobemer_output/ --k 15 \
           --strobe_w_min_offset 20 --strobe_w_max_offset 70 \
           --randstrobe_index
```

Minstrobes has the same parameters as hybridstrobes and randstrobes but are invoked with parameter `--minstrobe_index`


## Output

The output is a file `matches.tsv` in the output folder. You can se a custom outfile name with the parameter `--prefix`.
Output format is a tab separated file on the same format as MUMmer, with identical fields except the last one which is approximate reference sequence match length instead of what MUMmer produce:

```
>query_accession
ref_id  ref_pos query_pos   match_length_on_reference
```

Small example output from aligning sirv reads to transcripts (from the commands above) which also highlights the stobemers strength compared to kmers. While kmers can give a more nuanced differentiation (compare read hits to `SIRV606` and `SIRV616`) both the sequences are good candidates for downstream processing. In this small example, the strobemers produce fewer hits/less output needed for post clustering of matches, e.g., for downstream clustering/alignment/mapping. Notice that randstobe hit positions are currently not deterministic due to hash seed is set at each new pyhon instantiation. I will fix the hash seed in future implementations.


**Randstrobes (2,15,20,70)**
```
>41:650|d00e6247-9de6-485c-9b44-806023c51f13
SIRV606 35      92      487
SIRV616 35      92      473
>56:954|a23755a1-d138-489e-8efb-f119e679daf4
SIRV509 3       3       515
SIRV509 520     529     214
SIRV509 762     767     121
>106:777|0f79c12f-efed-4548-8fcc-49657f97a126
SIRV404 53      131     535
```

**kmers (k=30)**
```
>41:650|d00e6247-9de6-485c-9b44-806023c51f13
SIRV606 33      90      46
SIRV606 92      150     125
SIRV606 219     275     81
SIRV606 349     408     70
SIRV606 420     479     47
SIRV606 481     540     42
SIRV616 33      90      46
SIRV616 92      150     125
SIRV616 219     275     81
SIRV616 349     408     60
SIRV616 409     482     44
SIRV616 467     540     42
>56:954|a23755a1-d138-489e-8efb-f119e679daf4
SIRV509 68      72      141
SIRV509 230     233     100
SIRV509 331     335     105
SIRV509 435     442     40
SIRV509 475     483     36
SIRV509 579     585     41
SIRV509 621     627     46
SIRV509 695     701     44
SIRV509 812     815     53
>106:777|0f79c12f-efed-4548-8fcc-49657f97a126
SIRV404 53      131     58
SIRV404 128     208     127
SIRV404 283     364     30
SIRV404 422     494     142
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
