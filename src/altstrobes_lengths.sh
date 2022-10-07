#!/usr/bin/env bash 

for run in {1..1000}; do
  python3 matching_analysis_simulated.py --altstrobes_generalized
done

