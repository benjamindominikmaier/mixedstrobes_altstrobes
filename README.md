Mixedstrobes/Altstrobes/Multistrobes
===========

This repository contains information, scripts, and sourcecode to run the analysis for various strobemers in Maier and Sahlin, 2023 [bioRxiv LINK](https://www.biorxiv.org/content/10.1101/2022.10.13.512198). For info about strobemers, see [the repo](https://github.com/ksahlin/strobemers).


### What is a mixedstrobe?

Mixedstrobes is a mixed seeding technique that samples either a $k$-mer or a strobemer at a specified fraction. They may be sampled with all three strobemers seeding techniques (minstrobes, hybridstrobes and randstrobes) as well as altstrobes. Analogous to strobemers, we parameterize mixedstrobes as $(n,\ell,w_{min}, w_{max}, q)$, where $n$ is the number of strobes $n$, $\ell$ is the strobe length, $w_{min}$ and $w_{max}$ the minimum and maximum downstream offset to last window, and $q$ the strobemer fraction. Details on construction are found in [bioRxiv LINK](https://www.biorxiv.org/content/10.1101/2022.10.13.512198)

### What is an altstrobe?

Altstrobes are modified randstrobes where the strobe length is alternating between shorter and longer strobes. Contrary to randstrobes, where all strobes are of equal size (e.g. $k/2$ for strobemers of order 2), altstrobes are build up of one short strobe $k_s$ one longer strobe $k_l$, with $|k_s| + |k_l| = k$. We parameterize altstrobes as $(n,|k_s|,|k_l|, w_{min}, w_{max})$. Details on constructions are found in [bioRxiv LINK](https://www.biorxiv.org/content/10.1101/2022.10.13.512198)

### What is a multistrobe?

Multistrobes are generalized altstrobes where seeds with multiple different strobe lengths are generated within the same seeding pass. Contrary to altstrobes, where all strobes are of size $|k_s|,|k_l|$, multistrobes are selected in a range of lengths ranging from $|k_s|$ to $|k_l|=k-|k_s|$. We parameterize multistrobes analogously to altstrobes as $(n,|k_s|,|k_l|, w_{min}, w_{max})$. Details on constructions are found in [bioRxiv LINK](https://www.biorxiv.org/content/10.1101/2022.10.13.512198)

### This repository

The repository consists of 

- functions to generate strobemers in C++
- functions to generate strobemers in Python
- a tool `StrobeMap` implemented in both C++ and Python
- a modified version of 'minimap2' for strobemers
- scripts used for the evaluations in the [paper LINK](https://www.biorxiv.org/content/10.1101/2022.10.13.512198)

# C++ functions

The C++ library `strobemers_cpp/index.[cpp/hpp]` contains functions for creating altstrobes, multistrobes, mixedstrobes, randstobes, hybridstrobes, and minstrobes. 

You can copy the `index.cpp` and `index.hpp` files in the `strobemers_cpp` folder in this repository if you want to use either altstrobes, multistrobes, mixedstrobes, randstobes, hybridstrobes, or minstrobes in your project. They have the same user interface as described [here](https://github.com/ksahlin/strobemers#c-functions).


# Python functions

The `indexing_Maier_altstrobes.py` module located in the `modules` folder contains functions for generating k-mers, spaced k-mers, minimizers, and strobemers (minstrobes, hybridstrobes, randstrobes, mixedstrobes, altstrobes, and multistrobes) of any order. They have the same user interface as described [here](https://github.com/ksahlin/strobemers#python-functions).  

# StrobeMap (C++)


The tool `StrobeMap_Maier` is a program which roughly has the same interface as `MUMmer`. `StrobeMap_Maier` takes a reference and query file in multi-fasta or fastq format. It produces NAMs (Non-overlapping Approximate Matches) between the queries and references and outputs them in a format simular to nucmer/MUMmer. See [supplementary material](https://genome.cshlp.org/content/suppl/2021/10/19/gr.275648.121.DC1/Supplemental_Methods.pdf) Section A in the paper for definition of NAMs.

## Installation

`StrobeMap` (currently still without altstrobes, mixedstrobes and multistrobes) is available through [bioconda](https://bioconda.github.io/recipes/strobemap/README.html#package-strobemap). You can also acquire precompiled binaries for Linux from [here](https://github.com/benjamindominikmaier/mixedstrobes_altstrobes/tree/main/src/cpp/binaries). For example, for linux, simply do

```
wget ./mixedstrobes_altstrobes/src/cpp/binaries/StrobeMap
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
	-c Choice of protocol to use; kmers, minstrobes, hybridstrobes, altstrobes, randstrobes, multistrobes [randstrobes]. 
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

The output is a file `matches.tsv` in the output folder. You can set a custom outfile name with the parameter `--prefix`.
Output format is a tab separated file on the same format as MUMmer, with identical fields except the last one which is approximate reference sequence match length instead of what MUMmer produce:

```
>query_accession
ref_id  ref_pos query_pos   match_length_on_reference
```


# Minimap2-Strobemers

We implemented subsampled randstrobes, mixedstrobes, altstrobe, an multistrobe seeds in minimap2. A precompiled binary for Linux can be directly downloaded from [here](https://github.com/benjamindominikmaier/mixedstrobes_altstrobes/tree/main/src/minimap-strobemer/minimap2).

If you want to compile from the source, you need to have a C compiler, GNU make and zlib development files installed. Then type make in the source code directory to compile. If you see compilation errors, try make sse2only=1 to disable SSE4 code, which will make minimap2-strobemer slightly slower.

All default minimap2 settings can be found in the original minimap implementation manual [here](https://github.com/lh3/minimap2). Our implementation adds the following settings:
```
 --mixedstrobes  sample mixedrandstrobes80 instead of k-mers
 --altstrobes    sample altstrobes (k,2k) instead of k-mers
 --randstrobes   sample randstrobes (k,k) instead of k-mers
 -- multistrobes sample multistrobes (k_s, 2k-k_s) instead of k-mers
    -i INT       minimal strobe offset [25].
    -j INT       maximial strobe offset [50].
    -b INT       minimal k-mer size (multistrobes) [5].
```

# Evaluation Scripts

Are found in the `src` folder. Run instructions provided below

## matching_analysis_simulated.py

The script is used to compute the matching metrics (fraction of matches, the sequence coverage, the match coverage and the expected island size as definied in Sahlin 2021) of simulated sequences with various mutation rates and all strobemer methods.

When running the script without specifying any parameters (as in Maier & Sahlin, 2023), 1000 random DNA sequences (`--nr_exp`) of length 10,000nt (`--nr_exp`) are generated and subsequently mutated with mutation frequencies of 0.01, 0.05 and 0.1 ($--mut_freqs$) and equal chance for insertions, deletions and substitutions (see `--experiment_type`) to obtain corresponding query sequences. In the default mode, the sequences are seeded with all available methods (`kmers", "spaced_kmers_dense", "spaced_kmers_sparse", "minstrobes", "randstrobes", "hybridstrobes", "altstrobes", "mixedminstrobes", "mixedrandstrobes", "mixedhybridstrobes", "mixedaltstrobes`) and strobemer settings (2,15,25,50) (see `--k_size, --w, --orders, --w_low, --w_high`). By default, mixedstrobes are analyzed with strobe fractions ranging from 0.1 to 0.9, which can be changed using `--strobe_fractions`.

However, it is also possible to specify only one seeding technique using the `--method` parameter or sample generalized_altstrobes using `--altstrobes_generalized`. When sampling generalized altstrobes, the matching analysis is performed for altstrobes of all combinations from (1,k-1) to (k/2,k/2), whereby a lower boundary can be specified using `--k_boundary`, which is recommended as altstrobes with `k_s < 5` cause uniqueness issues and thus bad performance. Additionally, we implemented the possibility to sample e.g. strobemer-strobemer combinations using the `--mixedstrobes` parameter, which allows to sample a mix of any two techniques with distributions given by `--strobe_fractions`.

### Usage 

```
usage: matching_analysis_simulated.py [-h] [--L L] [--nr_exp NR_EXP]
                                      [--experiment_type EXPERIMENT_TYPE] [--subs_freq SUBS_FREQ]
                                      [--mut_freqs MUT_FREQS [MUT_FREQS ...]] [--k_size K_SIZE]
                                      [--w W] [--orders ORDERS [ORDERS ...]] [--w_low W_LOW]
                                      [--w_high W_HIGH]
                                      [--strobe_fractions STROBE_FRACTIONS [STROBE_FRACTIONS ...]]
                                      [--all_methods] [--method METHOD] [--altstrobes_generalized]
                                      [--altstrobes_size_distribution] [--multistrobes]
                                      [--multistrobes_size_distribution] [--k_boundary K_BOUNDARY]
                                      [--mixedstrobes]
                                      [--mixedstrobes_methods MIXEDSTROBES_METHODS] [--verbose]

Calc identity

optional arguments:
  -h, --help            show this help message and exit
  --L L                 Length of simulated sequences (default: 10000)
  --nr_exp NR_EXP       Number of simulated experiments (default: 1000)
  --experiment_type EXPERIMENT_TYPE
                        experiment type choose between "all", "controlled", "specified" or "only_subs" (default: all)
  --subs_freq SUBS_FREQ
                        substitution frequency among all mutations for --experiment_type "specified"; rest split evenly in insertions and deletions (default: 0.33)
  --mut_freqs MUT_FREQS [MUT_FREQS ...]
                        mutation frequencies [0,1] (default: [0.01, 0.05, 0.1])
  --k_size K_SIZE       k-mer/strobemer length (default: 30)
  --w W                 number of hashes used in a sliding window for thinning (w=1 means no thinning) (default: 1)
  --orders ORDERS [ORDERS ...]
                        List with orders of strobes to be analzyed (default: [2])
  --w_low W_LOW         minimum window offset to the previous window (wMin > 0) (default: 25)
  --w_high W_HIGH       maximum window offset to the previous window (wMin <= wMax) (default: 50)
  --strobe_fractions STROBE_FRACTIONS [STROBE_FRACTIONS ...]
                        Fraction of sampled strobemers, rest kmers (default: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
  --all_methods         perform matching analysis on simulated data for all (mixed-)strobemer and (mixed-)altstrobe seeding techniques (default: False)
  --method METHOD       choose seeding technique (default: none)
  --altstrobes_generalized
                        perform matching analysis on simulated data for altstrobes of all combinations from (1,k-1) to (k/2,k/2) (default: False)
  --altstrobes_size_distribution
                        perform matching analysis on simulated data for altstrobes with strobe size distribution (k_s, k_l) determined by strobe_fraction (default: False)
  --multistrobes        perform matching analysis on simulated data for multistrobes of all combinations from (k_boundary,k-k_boundary) to (k/2,k/2) (default: False)
  --multistrobes_size_distribution
                        perform matching analysis on simulated data for multistrobes of all combinations from (k_boundary,k-k_boundary) to (k/2,k/2) with strobe size distribution (k_s, k_l) determined
                        by strobe_fraction (default: False)
  --k_boundary K_BOUNDARY
                        minimum strobe length (k >= 4 recommended to ensure uniqueness) (default: 5)
  --mixedstrobes        perform matching analysis on simulated data for user defined mixed seeding techniques (default: False)
  --mixedstrobes_methods MIXEDSTROBES_METHODS
                        List with two seeding methods to sample mixedstrobes (default: ['randstrobes', 'kmers'])
  --verbose
```

### Output

The output format is as follows:
```
Method & Setting & Fraction of Matches & Sequence Coverage & Match Coverage & Expected Island Size & Mutation Frequency, e.g.
minstrobes  &  (2, 15, 25, 50)  &  67.0 & 75.7 & 96.2 & 2.1  &  0.01
```

## matching_analysis_bio

The script is used to compute the matching metrics (as definied in Sahlin 2021 [LINK](https://www.genome.org/cgi/doi/10.1101/gr.275648.121)) of biological sequences for all strobemer methods. The query sequences (`--queries`) are split up in 
disjoint segments of length (`--segment`) and mapped to the reference (`--references`) before the collinear chain solution of raw unmerged hits is determined for each segment and the matching metrics computed from it. The collinear chain solution takes only the longest collinear chain of hits into account, thus assuming the most likely location and avoiding to overcount "spurious" hits (see Sahlin, 2021 [LINK](https://www.genome.org/cgi/doi/10.1101/gr.275648.121)). Details on implementation are found in [bioRxiv LINK](https://www.biorxiv.org/content/10.1101/2022.10.13.512198) Supplementary Section S4.

### Usage
```
usage: matching_analysis_bio.py [-h] [--queries QUERIES] [--references REFERENCES] [--k K] [--strobe_w_min_offset STROBE_W_MIN_OFFSET] [--strobe_w_max_offset STROBE_W_MAX_OFFSET] [--w W] [--n N] [--strobe_fraction STROBE_FRACTION] [--dont_merge_matches] [--outfolder OUTFOLDER] [--prefix PREFIX] [--kmer_index] [--minstrobe_index] [--randstrobe_index] [--hybridstrobe_index] [--altstrobe_index] [--multistrobe_index] [--mixedminstrobe_index] [--mixedrandstrobe_index] [--mixedhybridstrobe_index] [--mixedaltstrobe_index] [--segment SEGMENT] [--selfalign] [--rev_comp]

Calc identity

optional arguments:
  -h, --help            show this help message and exit
  --queries QUERIES     Path to query fasta or fastq file (default: False)
  --references REFERENCES
                        Path to reference fasta or fastq file (default: False)
  --k K                 Length of kmer/all strobes combined (default: 15)
  --strobe_w_min_offset STROBE_W_MIN_OFFSET
                        Strobemer window start offset from first k-mer. If kmer start at pos i, first window will start at i+strobe_w_min_offset. Default: 20nt donwstream from start of first kmer.
                        (default: 20)
  --strobe_w_max_offset STROBE_W_MAX_OFFSET
                        Strobemer window end. If kmer start at pos i, first window will stop at i+strobe_w_max_offset. Default: 70nt donwstream from start of first kmer. (default: 70)
  --w W                 Thinning window size applied to reference sequences (default = 1, i.e., no thinning) (default: 1)
  --n N                 Order on strobes (default: 2)
  --strobe_fraction STROBE_FRACTION
                        Fraction of sampled strobemers, rest kmers (default: 1)
  --outfolder OUTFOLDER
                        Folder to output TSV match file. (default: output_matching_analysis_bio)
  --prefix PREFIX       Filename prefix (default "matches"). (default: matches)
  --kmer_index          Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for k-mers (default: False)
  --minstrobe_index     Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for minstrobes (default: False)
  --randstrobe_index    Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for randstrobes (default: False)
  --hybridstrobe_index  Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for hybridstrobes (default: False)
  --altstrobe_index     Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for altstrobes (default: False)
  --multistrobe_index   Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for multistrobes (default: False)
  --mixedminstrobe_index
                        Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for mixed minstrobes/kmers based on
                        --strobe_fraction (default: False)
  --mixedrandstrobe_index
                        Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for mixed randstrobes/kmers based on
                        --strobe_fraction (default: False)
  --mixedhybridstrobe_index
                        Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for mixed hybridstrobes/kmers based on
                        --strobe_fraction (default: False)
  --mixedaltstrobe_index
                        Produce chains of matches that are in identical order in both sequences (collinear chaining algorithm) and compute matching metrics for mixed altstrobes/kmers based on
                        --strobe_fraction (default: False)
  --segment SEGMENT     segment length for computing the collinear chain solution of the raw hits (default: 2000)
  --selfalign           Aligns sequences to itself (mainly used for bugfixing). Default is not align sequences to themselves if the same file is given as references and queries. (default: False)
  --rev_comp            Match reverse complement of reads (output to separate file) (default: False)
```
### Output

The script creates two output files, whereby one contains information about all segments from each query:
```
altstrobe.txt

Query Number & Query Length & Segment Number & Segment Length & #Matches & Fraction of Matches & Sequence Coverage & Match Coverage & Expected Island Size & Collinear E-Size
1 & 52197 & 1 & 2000 & 264 & 0.03 & 0.41 & 0.55 & 53.99 & 25.1
1 & 52197 & 2 & 2000 & 318 & 0.02 & 0.43 & 0.6 & 52.98 & 21.0
```

and the second summary file contains the average matching metrics (weighted mean over all segments by segment length) for each of the queries:


```
altstrobe (summary).txt

Query Number & Query Length & Fraction of Matches & Sequence Coverage & Match Coverage & Expected Island Size
>Query 1 & 52197 & 2.5 & 38.4 & 52.5 & 147.1 & 555.8
>Query 2 & 51992 & 3.4 & 51.4 & 70.3 & 35.8 & 747.5
```

## uniqueness_analysis

This script computes the uniqueness of seeds and the expected number of hits (E-size) for either a given reference sequence or a simulated sequence. E-size is a measure of how repetitive the seeds in a query sequence are, on average, in a reference dataset; details on calculation are given in Sahlin, 2022 [LINK](https://doi-org.ezp.lib.cam.ac.uk/10.1186/s13059-022-02831-7).

### Usage

```
usage: uniqueness_analysis.py [-h] [--fasta FASTA] [--kmers] [--minstrobes] [--mixedminstrobes] [--randstrobes] [--mixedrandstrobes] [--hybridstrobes] [--mixedhybridstrobes] [--spaced_dense]
                              [--spaced_sparse] [--altstrobes] [--multistrobes] [--mixedaltstrobes] [--order ORDER] [--strobe_fraction STROBE_FRACTION] [--k_sizes K_SIZES [K_SIZES ...]]
                              [--k_boundary K_BOUNDARY] [--w W] [--w_low W_LOW] [--w_high W_HIGH] [--altstrobes_generalized ALTSTROBES_GENERALIZED] [--L L] [--nr_exp NR_EXP]

Calc identity

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA         Path to consensus fastq file(s) (default: False)
  --kmers               Seeding k-mers (default: False)
  --minstrobes          Seeding minstrobes (default: False)
  --mixedminstrobes     Seeding mixedminstrobes (minstrobes/k-mers) with a user-defined --strobe_fraction (default: False)
  --randstrobes         Seeding randstrobes (default: False)
  --mixedrandstrobes    Seeding mixedrandstrobes (randstrobes/k-mers) with a user-defined --strobe_fraction (default: False)
  --hybridstrobes       Seeding hybridstrobes (default: False)
  --mixedhybridstrobes  Seeding mixedhybridstrobes (hybridstrobes/k-mers) with a user-defined --strobe_fraction (default: False)
  --spaced_dense        Seeding spaced k-mers (dense) (default: False)
  --spaced_sparse       Seeding spaced k-mers (sparse) (default: False)
  --altstrobes          Seeding altstrobes (default: False)
  --multistrobes        Seeding multistrobes of all combinations from (k_boundary,k-k_boundary) to (k/2,k/2) (default: False)
  --mixedaltstrobes     Seeding mixedaltstrobes (altstrobes/k-mers) with a user-defined --strobe_fraction (default: False)
  --order ORDER         Order on strobes (default: 2)
  --strobe_fraction STROBE_FRACTION
                        Fraction of sampled strobemers, rest kmers (default: 1)
  --k_sizes K_SIZES [K_SIZES ...], --nargs-int-type K_SIZES [K_SIZES ...]
                        List with strobe lengths to be analyzed (default: [18, 24, 30, 36])
  --k_boundary K_BOUNDARY
                        minimum strobe length (k >= 4 recommended to ensure uniqueness) (default: 5)
  --w W                 number of hashes used in a sliding window for thinning (w=1 means no thinning) (default: 1)
  --w_low W_LOW         minimum window offset to the previous window (wMin > 0) (default: 25)
  --w_high W_HIGH       maximum window offset to the previous window (wMin <= wMax) (default: 50)
  --altstrobes_generalized ALTSTROBES_GENERALIZED
                        Choose k-size to seed altstrobes with strobe combinations from (1,k-1) to (k-1,1) (default: 0)
  --L L                 Length of simulated sequences (only if no fasta file is provided) (default: 100000)
  --nr_exp NR_EXP       Number of simulated experiments (only if no fasta file is provided) (default: 1)

```

### Output

The output format is as follows:
seeding_technique, k_size, accession_number (only for biological sequences), mean, median, lower_75, upper_75, percent_unique, ehits

```
python3 uniqueness_analysis.py --randstrobes

Simulating sequence of length:  10000

randstrobes2,18,,1.0,1,1,1,100.0,1.0
randstrobes2,24,,1.0,1,1,1,100.0,1.0
randstrobes2,30,,1.0,1,1,1,100.0,1.0
randstrobes2,36,,1.0,1,1,1,100.0,1.0

```

## generate_seed

This script samples $n$ reads (`nr_exp`) from a reference (`--reference`) and mutates them for the different experimental conditions (`--mut_freqs`). The sampled reads (unmutated and mutated) are saved to different output files (`--prefix`) in the outfolder (`--outfolder`). The reads are sampled from random positions of the reference and it is possible to specify whether sampled reads may contain unspecified nucleotides using `--tolerate_N`. In our analysis, this script is used to produce simulated reads for the minimap analysis (see next section).

### Usage

```
usage: generate_seeds.py [-h] [--L L] [--nr_exp NR_EXP] [--experiment_type EXPERIMENT_TYPE] [--mut_freqs MUT_FREQS] [--reference REFERENCE] [--rev_comp] [--outfolder OUTFOLDER] [--prefix PREFIX] [--tolerate_N]

Calc identity

optional arguments:
  -h, --help            show this help message and exit
  --L L                 Length of simulated sequences (default: 10000)
  --nr_exp NR_EXP       Number of simulated experiments (default: 100000)
  --experiment_type EXPERIMENT_TYPE
                        experiment type choose between "all", "controlled or "only_subs" (default: all)
  --mut_freqs MUT_FREQS
                        mutation frequencies [0.01, 0.05, 0.1] (default: [0.01, 0.05, 0.1])
  --reference REFERENCE
                        path to fasta file with reference genome to generate reads from (default: None)
  --rev_comp            Sampling reads from both strands (default: False)
  --outfolder OUTFOLDER
                        Folder to output FA read file. (default: ../output/simulated_reads)
  --prefix PREFIX       Filename prefix (default "matches"). (default: reads)
  --tolerate_N          Allow Ns in simulated reads (default: False)
```

### Output

The output is a file `reads.txt` containing the unmutated sampled reads and multiple `reads_<mutation_frequency>_mutated.txt` files with the corresponding mutated reads in the output folder. You can set a custom outfile name with the parameter `--prefix`. Output format is fasta.

## analyze_minimap

This script takes a SAM output file from minimap2-strobemers (`--input_file`) and computes the fraction of correctly mapped reads (as defined by the existence of an overlap between query and reference position and right direction) and the average (mean) number of correctly mapped nucleotides. 

### Usage

```
usage: analyze_minimap.py [-h] [--input_file INPUT_FILE] [--L L] [--N N]

Calc identity

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Name of SAM output file from minimap2 (default: None)
  --L L                 Length of each query read (default: 10000)
  --N N                 Number of query reads (default: 100000)
```

### Output

This script outputs the fraction of correctly mapped reads (1st line) and the average (mean) number of correctly mapped nucleotides (2nd line) as well as the number of incorrectly mapped reads on the forward strand, the number of incorrectly mapped reads on the reverse strand and the total number of mapped reads (3rd line) to the command line.

```
python3 analyze_minimap.py --input_file ../output/minimap/cmh13_kmers15_00.sam --L 10000 --N 100000

0.99513
9950.88829
238 249 100000
```

## Stochasticity in Seed Construct

This script performs simulations showing how stochasticity in seed construct influence probability of `w` consecutive seeds producing at least one match in a region of length `2w = 128` between sequences for various seed constructs.

### Usage (fig2A_analysis.py)

```
usage: fig2A_analysis.py [-h] [--N_SIM N_SIM] [--thinning THINNING] k w_min w_max outcsv outplot_prefix

Calc identity

positional arguments:
  k               k-mer size
  w_min           Window size
  w_max           Window size
  outcsv          output CSV file for plotting.
  outplot_prefix  output pdf file prefix.

optional arguments:
  -h, --help      show this help message and exit
  --N_SIM N_SIM   Number of simulations for p_matches experiment (default: 1000)
  --thinning      Thinning level, 1/X. (default X=1 no thinning)
```

### Output

This script produces csv.files with the following columns: `w_min,w_max,seeding,N_m,stddev_est,m` as well as corresponding plots.

### Usage (fig2B-C_conditional_entropy_on_c.py)

```
usage: fig2B-C_conditional_entropy_on_c.py [-h] [--max_mut_freq MAX_MUT_FREQ] [--N_SIM N_SIM] k w_min w_max outcsv outplot_prefix


Calc identity

positional arguments:
  k               k-mer size
  w_min           Window size
  w_max           Window size
  outcsv          output CSV file for plotting.
  outplot_prefix  output pdf file prefix.

optional arguments:
  -h, --help      show this help message and exit
  --N_SIM N_SIM   Number of simulations for p_matches experiment (default: 1000)
  --max_mut_freq  Maximum mutation frequency for simulations
```

### Output

This script produces csv.files with the following columns: `w_min,w_max,type,x,m,y,mean_error_y_est,analysis,ymetric` as well as corresponding plots. The y-metric is either `p_match` or `Entropy`.

## Figures in Maier & Sahlin, 2023

All figures in Maier & Sahlin, 2023 can be generated running the R-notebook `figures.rmd`. All data tables as well as all figures used in this manuscript can be found in the output folder in this repository.

CREDITS
----------------

Benjamin D. Maier and Kristoffer Sahlin, Entropy predicts sensitivity of pseudo-random seeds, bioRxiv, January 2023; https://www.biorxiv.org/content/10.1101/2022.10.13.512198
