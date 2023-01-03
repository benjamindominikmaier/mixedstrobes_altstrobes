#!/usr/bin/env bash 

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --kmers
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --spaced_dense
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --spaced_sparse

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --minstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --randstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --hybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --minstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --randstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --hybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --altstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --multistrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.1 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.1 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.1 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.1 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.1 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.1 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.1 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.2 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.2 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.2 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.2 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.2 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.2 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.2 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.3 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.3 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.3 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.3 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.3 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.3 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.3 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.4 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.4 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.4 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.4 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.4 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.4 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.4 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.5 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.5 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.5 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.5 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.5 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.5 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.5 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.6 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.6 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.6 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.6 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.6 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.6 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.6 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.7 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.7 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.7 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.7 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.7 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.7 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.7 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.8 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.8 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.8 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.8 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.8 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.8 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.8 --mixedhybridstrobes

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.9 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.9 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.9 --mixedhybridstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.9 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.9 --mixedminstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.9 --mixedrandstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 3 --strobe_fraction 0.9 --mixedhybridstrobes
