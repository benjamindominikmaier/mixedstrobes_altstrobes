#!/usr/bin/env bash 

python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.1 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.3 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.5 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.7 --mixedaltstrobes
python3 uniqueness_analysis.py --fasta ../data/chr21.fna --order 2 --strobe_fraction 0.9 --mixedaltstrobes

