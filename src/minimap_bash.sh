#!/usr/bin/env bash 

#kmers
./minimap2 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads > cmh13_kmers15_00.sam
./minimap2 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.01_mutated > cmh13_kmers15_01.sam
./minimap2 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.05_mutated > cmh13_kmers15_05.sam
./minimap2 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.1_mutated > cmh13_kmers15_10.sam
./minimap2 -k 28 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads > cmh13_kmers28_00.sam
./minimap2 -k 28 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.01_mutated > cmh13_kmers28_01.sam
./minimap2 -k 28 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.05_mutated > cmh13_kmers28_05.sam
./minimap2 -k 28 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.1_mutated > cmh13_kmers28_10.sam

#altstrobes
./minimap2 --altstrobes -k 9 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads > cmh13_altstrobes9_00.sam
./minimap2 --altstrobes -k 9 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.01_mutated > cmh13_altstrobes9_01.sam
./minimap2 --altstrobes -k 9 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.05_mutated > cmh13_altstrobes9_05.sam
./minimap2 --altstrobes -k 9 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.1_mutated > cmh13_altstrobes9_10.sam

#randstrobes
./minimap2 --randstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads > cmh13_randstrobes14_00.sam
./minimap2 --randstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.01_mutated > cmh13_randstrobes14_01.sam
./minimap2 --randstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.05_mutated > cmh13_randstrobes14_05.sam
./minimap2 --randstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.1_mutated > cmh13_randstrobes14_10.sam

#mixedstrobes
./minimap2 --mixedstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads > cmh13_mixedstrobes14_00.sam
./minimap2 --mixedstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.01_mutated > cmh13_mixedstrobes14_01.sam
./minimap2 --mixedstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.05_mutated > cmh13_mixedstrobes14_05.sam
./minimap2 --mixedstrobes -k 14 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.1_mutated > cmh13_mixedstrobes14_10.sam

#multistrobes
./minimap2 --multistrobes -k 14 -b 6 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads > cmh13_multistrobes9_00.sam
./minimap2 --multistrobes -k 14 -b 6 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.01_mutated > cmh13_multistrobes9_01.sam
./minimap2 --multistrobes -k 14 -b 6 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.05_mutated > cmh13_multistrobes9_05.sam
./minimap2 --multistrobes -k 14 -b 6 -a ../data/chm13v2.0.fa ../output/simulated_reads/reads_0.1_mutated > cmh13_multistrobes9_10.sam
