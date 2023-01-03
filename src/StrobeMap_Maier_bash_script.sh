#!/usr/bin/env bash 

# kmers
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 30 --kmer_index --rev_comp --prefix kmer --outfolder ../output/data\ tables/StrobeMap/


# stromebers
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --randstrobe_index --rev_comp --prefix randstrobe2 --outfolder ../output/data\ tables/StrobeMap/

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --randstrobe_index --rev_comp --prefix randstrobe3 --outfolder ../output/data\ tables/StrobeMap/


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --minstrobe_index --rev_comp --prefix minstrobe2 --outfolder ../output/data\ tables/StrobeMap/

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --minstrobe_index --rev_comp --prefix minstrobe3 --outfolder ../output/data\ tables/StrobeMap/


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --hybridstrobe_index --rev_comp --prefix hybridstrobe2 --outfolder ../output/data\ tables/StrobeMap/

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --hybridstrobe_index --rev_comp --prefix hybridstrobe3 --outfolder ../output/data\ tables/StrobeMap/


# mixedstrobemers
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.1 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.1
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.1 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.1

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.1 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.1
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.1 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.1

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.1 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.1
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.1 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.1


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.2 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.2
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.2 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.2

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.2 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.2
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.2 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.2

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.2 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.2
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.2 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.2


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.3 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.3
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.3 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.3

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.3 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.3
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.3 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.3

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.3 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.3
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.3 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.3


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.4 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.4
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.4 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.4

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.4 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.4
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.4 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.4

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.4 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.4
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.4 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.4


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.5 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.5
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.5 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.5

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.5 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.5
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.5 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.5

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.5 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.5
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.5 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.5


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.6 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.6
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.6 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.6

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.6 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.6
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.6 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.6

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.6 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.6
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.6 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.6


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.7 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.7
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.7 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.7

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.7 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.7
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.7 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.7

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.7 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.7
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.7 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.7


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.8 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.8
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.8 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.8

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.8 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.8
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.8 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.8

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.8 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.8
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.8 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.8


python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe2_0.9 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.9
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedrandstrobe_index --rev_comp --prefix mixedrandstrobe3_0.9 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.9

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe2_0.9 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.9
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedminstrobe_index --rev_comp --prefix mixedminstrobe3_0.9 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.9

python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 15 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe2_0.9 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.9
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 3 --strobe_w_min_offset 11 --strobe_w_max_offset 100 --mixedhybridstrobe_index --rev_comp --prefix mixedhybridstrobe3_0.9 --outfolder ../output/data\ tables/StrobeMap/
 --strobe_fraction 0.9
 
# multistrobes
python3 StrobeMap_Maier --queries ../data/GCF_003018135.1_ASM301813v1_genomic.fna --references ../data/GCF_003018575.1_ASM301857v1_genomic.fna --w 1 --k 10 --n 2 --strobe_w_min_offset 16 --strobe_w_max_offset 100 --multistrobe_index --rev_comp --prefix multistrobes2 --outfolder ../output/data\ tables/StrobeMap/
