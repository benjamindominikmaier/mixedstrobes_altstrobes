#!/usr/bin/env python3.9
# -*- coding: utf-8 -*

"""Mash distance calculation of simulated sequences for ANI estimation.

Please note, that the ANI estimations should be adjusted to spaced kmers and
strobemers by applying the following correction:
ANI = D + cf * (1-D)
with correction factor cf = 0.075 (kmers), 0.165 (strobemers), 0.3/0.6 (spaced kmers)

Approach based on Ondov, B.D., Treangen, T.J., Melsted, P. et al.
"Mash: fast genome and metagenome distance estimation using MinHash."
Genome Biol 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x
"""

__authors__ = ["Benjamin D. Maier"]
__copyright__ = "Copyright Benjamin D. Maier & Kristoffer Sahlin | Sahlin Group"
__organization__ = "Department of Mathematics, Science for Life Laboratory, Stockholm University, 106 91, Stockholm, Sweden."
__credits__ = ["Benjamin D. Maier & Kristoffer Sahlin"]
__contact__ = "bmaier [at] ebi.ac.uk"
__date__ = "2023/03/10"
__created__ = "2023/02/27"
__deprecated__ = False
__license__ = "MIT"
__maintainer__ = "Kristoffer Sahlin"
__email__ = "kristoffer.sahlin [at] scilifelab.se"
__status__ = "DSML Lvl. 1 - Concept"

import os
import sys
import argparse
import random
import numpy as np
import collections
import math
import statistics
from typing import Iterator
from itertools import chain

from modules import help_functions
from modules import indexing_Maier_altstrobes as indexing


def analyze_strobemers(seq1: str, seq2: str, k_size: int, order: int,
                       hash_fcn: str, w: int, w_low: int = 0, w_high: int = 50,
                       fraction: float = 0.5) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), all_pos_vectors (list),
              match coverage (float)
    """
    mixed_methods = (
        "mixedminstrobes", "mixedminstrobes_generalized", "mixedrandstrobes",
        "mixedrandstrobes_generalized", "mixedhybridstrobes",
        "mixedhybridstrobes_generalized"
    )

    assert k_size % order == 0, "Not even kmer length, results will be different"
    if hash_fcn in mixed_methods:
        strobemers1 = getattr(indexing, hash_fcn)(seq1, k_size, w_low, w_high, w, order, fraction)
        strobemers2 = getattr(indexing, hash_fcn)(seq2, k_size, w_low, w_high, w, order, fraction)
    else:
        strobemers1 = getattr(indexing, hash_fcn)(seq1, k_size, w_low, w_high, w, order)
        strobemers2 = getattr(indexing, hash_fcn)(seq2, k_size, w_low, w_high, w, order)

    matches = len(set(strobemers1.values()) & set(strobemers2.values()))
    non_matches = len(strobemers1) + len(strobemers2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI


def analyze_altstrobes_generalized(seq1: str, seq2: str, k_size1: int,
                                   k_size2: int, order: int, w: int,
                                   w_low: int = 0, w_high: int = 50,
                                   fraction: float = 0.5) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated altstrobe seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), match coverage (float)
    """
    altstrobes1 = indexing.altstrobes_generalized(seq1, k_size1, k_size2, w_low, w_high, w, order, arg=0)
    altstrobes2 = indexing.altstrobes_generalized(seq2, k_size1, k_size2, w_low, w_high, w, order, arg=0)
    altstrobes2_rev = indexing.altstrobes_generalized(seq2, k_size1, k_size2, w_low, w_high, w, order, arg=1)

    matches = len(set(altstrobes1.values()) & (set(altstrobes2.values()) | set(altstrobes2_rev.values())))
    non_matches = len(altstrobes1) + len(altstrobes2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI


def analyze_multistrobes(seq1: str, seq2: str, k_size: int, order: int, w: int,
                         w_low: int = 0, w_high: int = 50, fraction: float = 0.5,
                         k_boundary: int=5, arg: int=50) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated altstrobe seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), match coverage (float)
    """
    multistrobes1 = indexing.multistrobes(seq1, k_size, w_low, w_high, w, order, k_boundary=k_boundary, arg=arg)
    multistrobes2 = indexing.multistrobes(seq2, k_size, w_low, w_high, w, order, k_boundary=k_boundary, arg=arg)

    matches = len(set(multistrobes1.values()) & (set(multistrobes2.values())))
    non_matches = len(multistrobes1) + len(multistrobes2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI


def analyze_kmers(seq1: str, seq2: str, k_size: int, w: int) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated k-mer seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of the kmers
    :param w: number of kmers used in a sliding window for thinning (w=1 means no thinning)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), all_pos_vectors (list),
              match coverage (float)
    """
    # kmers
    kmers_pos1 = indexing.kmers(seq1, k_size, w)
    kmers_pos2 = indexing.kmers(seq2, k_size, w)

    matches = len(set(kmers_pos1.values()) & set(kmers_pos2.values()))
    non_matches = len(kmers_pos1) + len(kmers_pos2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI


def analyze_spaced_kmers(seq1: str, seq2: str, k_size: int, span_size: int, w: int) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated spaced kmer seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of the kmers
    :param span_size: length between first and last position
    :param w: number of spaced_kmers used in a sliding window for thinning (w=1 means no thinning)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), all_pos_vectors (list),
              match coverage (float)
    """
    positions = set(random.sample(range(1, span_size-1), k_size-2))
    positions.add(0)
    positions.add(span_size - 1)  # asserts first and last position is sampled so that we have a spaced kmer of length span size
    spaced_kmers_seq1 = indexing.spaced_kmers(seq1, k_size, span_size, positions, w)
    spaced_kmers_seq2 = indexing.spaced_kmers(seq2, k_size, span_size, positions, w)
    matches = len(set(spaced_kmers_seq1.values()) & set(spaced_kmers_seq2.values()))
    non_matches = len(spaced_kmers_seq1) + len(spaced_kmers_seq2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI

def analyze_altstrobes(seq1: str, seq2: str, k_size: int, order: int,
                       hash_fcn: str, w: int, w_low: int = 0, w_high: int = 50,
                       fraction: float = 0.5, arg: int = 50) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated altstrobe seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), all_pos_vectors (list),
              match coverage (float)
    """

    assert k_size % int(1.5*order) == 0, "Not even kmer length, results will be different"
    altstrobes1 = indexing.altstrobes(seq1, k_size, w_low, w_high, w, order, arg=arg)
    altstrobes2 = indexing.altstrobes(seq2, k_size, w_low, w_high, w, order, arg=arg)

    order = 1.5*order
    matches = len(set(altstrobes1.values()) & set(altstrobes2.values()))
    non_matches = len(altstrobes1) + len(altstrobes2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI


def analyze_mixedaltstrobes(seq1: str, seq2: str, k_size: int, order: int,
                            hash_fcn: str, w: int, w_low: int = 0,
                            w_high: int = 50, fraction: float = 0.5) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated mixedaltstrobe seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence whichfrom itertools import chain is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), all_pos_vectors (list),
              match coverage (float)
    """

    assert k_size % int(1.5*order) == 0, "Not even kmer length, results will be different"
    altstrobes1 = indexing.mixedaltstrobes(seq1, k_size, w_low, w_high, w, order, fraction)
    altstrobes2 = indexing.mixedaltstrobes(seq2, k_size, w_low, w_high, w, order, fraction)

    order = 1.5*order  # increase order by 1 as information about short-long combination requires one extra order

    matches = len(set(altstrobes1.values()) & set(altstrobes2.values()))
    non_matches = len(altstrobes1) + len(altstrobes2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI


def analyze_mixedstrobes(method1: str, method2:str, seq1: str, seq2: str,
                         k_size: int, order: int, hash_fcn: str, w: int,
                         w_low: int = 0, w_high: int = 50,
                         fraction: float = 0.5) -> tuple:
    """
    Computes number of matches, fraction of matches, sequence coverage,
    match coverage and expected island size for generated mixed strobemer seeds

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence whichfrom itertools import chain is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: number of matches (int), fraction of matches (float),
              sequence coverage (float), gap lengths (list), all_pos_vectors (list),
              match coverage (float)
    """

    assert k_size % order == 0, "Not even kmer length, results will be different"
    mixedstrobes1 = indexing.mixedstrobes(method1, method2, seq1, k_size, w_low, w_high, w, order, fraction)
    mixedstrobes2 = indexing.mixedstrobes(method1, method2, seq2, k_size, w_low, w_high, w, order, fraction)

    matches = len(set(mixedstrobes1.values()) & set(mixedstrobes2.values()))
    non_matches = len(mixedstrobes1) + len(mixedstrobes2) - matches
    if matches == 0:
        return 0
    else:
        ANI = 1 + 1/k_size * math.log((2*(matches/non_matches))/(1+(matches/non_matches)))
        return ANI


def main(args):
    """
    Computes matching metrics (m, mp, sc, mc, e-size) for simulated sequences
    using various seeding techniques (kmer-based, strobemers, altstrobes)
    """
    log_file = open('../output/data tables/ANI_logger','a+')

    for mut_freq in args.mut_freqs:
        if args.experiment_type == "specified":
            print("EXPERIMENT TYPE: specified ({0}% subs)".format(100*args.subs_freq))
        print("MUTATION RATE:", mut_freq)
        results = dict()

        for exp_id in range(args.nr_exp):
            if ((args.verbose) & (exp_id % 100 == 0)):
                print("Analyzed {0} simulated experiments". format(exp_id))

            seq1 = "".join([random.choice("ACGT") for i in range(args.L)])

            # controlled or random experiment
            if args.experiment_type == 'only_subs':
                muts = set(random.sample(range(len(seq1)), int(args.L*mut_freq)))
                seq2 = "".join([
                    seq1[i] if i not in muts
                    else random.choice([help_functions.reverse_complement(seq1[i])])
                    for i in range(len(seq1))
                ])
            elif args.experiment_type == 'controlled':
                # muts = set(range(15,L,15)) # every 15th nt for figure 2 only!
                muts = set(range(20, args.L, 20))
                seq2 = "".join([
                    seq1[i] if i not in muts
                    else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")])
                    for i in range(len(seq1))
                ])
            elif args.experiment_type == 'all':
                muts = set(random.sample(range(len(seq1)), int(args.L*mut_freq)))
                seq2 = "".join([
                    seq1[i] if i not in muts
                    else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")])
                    for i in range(len(seq1))
                ])
            elif args.experiment_type == 'specified':
                muts = set(random.sample(range(len(seq1)), int(args.L*mut_freq)))
                subs = set(random.sample(muts, int(len(muts)*args.subs_freq)))
                seq2 = "".join([
                    seq1[i] if i not in muts
                    else random.choice([help_functions.reverse_complement(seq1[i])]) if i in subs
                    else random.choice(['', seq1[i] + random.choice("ACGT")])
                    for i in range(len(seq1))
                ])

            else:
                print("Wrong experiment label specified")
                sys.exit()

            for hash_fcn in args.methods:
                if hash_fcn == "kmers":
                    results.setdefault("kmers", {"ANI": []})
                    ANI = analyze_kmers(seq1, seq2, args.k_size, args.w)
                    results["kmers"]["ANI"].append(ANI)

                elif hash_fcn == "spaced_kmers_dense":
                    results.setdefault("spaced_kmers_dense", {"ANI": []})
                    ANI = analyze_spaced_kmers(seq1, seq2, args.k_size, args.k_size+args.k_size//2, args.w)
                    results["spaced_kmers_dense"]["ANI"].append(ANI)

                elif hash_fcn == "spaced_kmers_sparse":
                    results.setdefault("spaced_kmers_sparse", {"ANI": []})
                    ANI = analyze_spaced_kmers(seq1, seq2, args.k_size, 3*args.k_size, args.w)
                    results["spaced_kmers_sparse"]["ANI"].append(ANI)

                elif hash_fcn == "altstrobes_size_distribution":
                    results.setdefault("altstrobes_size_distribution", dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results["altstrobes_size_distribution"].setdefault((order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction), {"ANI": []})
                            ANI = analyze_altstrobes(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, arg=100*fraction)
                            results["altstrobes_size_distribution"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["ANI"].append(ANI)

                elif hash_fcn == "altstrobes":
                    results.setdefault("altstrobes", dict())
                    for order in args.orders:
                        results["altstrobes"].setdefault((order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high), {"ANI": []})
                        ANI = analyze_altstrobes(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high)
                        results["altstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high)]["ANI"].append(ANI)

                elif hash_fcn == "altstrobes_generalized":
                    results.setdefault("altstrobes_generalized", dict())
                    for k1, k2 in args.strobe_lengths:
                        results["altstrobes_generalized"].setdefault((k1, k2), {"ANI": []})
                        ANI = analyze_altstrobes_generalized(seq1, seq2, k1, k2, 2, args.w, w_low=args.w_low, w_high=args.w_high)
                        results["altstrobes_generalized"][(k1, k2)]["ANI"].append(ANI)

                elif hash_fcn == "mixedaltstrobes":
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results["mixedaltstrobes"].setdefault((order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction), {"ANI": []})
                            ANI = analyze_mixedaltstrobes(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, fraction=fraction)
                            results["mixedaltstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["ANI"].append(ANI)

                elif hash_fcn == "multistrobes":
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        results["multistrobes"].setdefault((order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1), {"ANI": []})
                        ANI = analyze_multistrobes(seq1, seq2, args.k_size, order, args.w, w_low=args.w_low, w_high=args.w_high, k_boundary=args.k_boundary)
                        results["multistrobes"][(order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1)]["ANI"].append(ANI)

                elif hash_fcn == "multistrobes_size_distribution":
                    results.setdefault("multistrobes_size_distribution", dict())
                    for fraction in args.strobe_fractions:
                        results["multistrobes_size_distribution"].setdefault((2, args.k_size, args.w_low, args.w_high, fraction), {"ANI": []})
                        ANI = analyze_multistrobes(seq1, seq2, args.k_size, 2, args.w, w_low=args.w_low, w_high=args.w_high, k_boundary=args.k_boundary, arg=100*fraction)
                        results["multistrobes_size_distribution"][(2, args.k_size, args.w_low, args.w_high, fraction)]["ANI"].append(ANI)

                elif hash_fcn == "mixedstrobes":
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results["mixedstrobes"].setdefault((args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction), {"ANI": []})
                            ANI = analyze_mixedstrobes(args.method1, args.method2, seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, fraction=fraction)
                            results["mixedstrobes"][(args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction)]["ANI"].append(ANI)

                elif hash_fcn in ("mixedminstrobes", "mixedminstrobes_generalized", "mixedrandstrobes", "mixedrandstrobes_generalized", "mixedhybridstrobes", "mixedhybridstrobes_generalized"):
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results[hash_fcn].setdefault((order, args.k_size//order, args.w_low, args.w_high, fraction), {"ANI": []})
                            ANI = analyze_strobemers(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, fraction=fraction)
                            results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high, fraction)]["ANI"].append(ANI)

                else:
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        results[hash_fcn].setdefault((order, args.k_size//order, args.w_low, args.w_high), {"ANI": []})
                        ANI = analyze_strobemers(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high)
                        results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high)]["ANI"].append(ANI)

        for protocol in results:
            if protocol == "kmers" or protocol == "spaced_kmers_sparse" or protocol == "spaced_kmers_dense":
                log_file.write(protocol + " & - & " + ",".join([str(round(r, 4)) for r in results[protocol]["ANI"]]) + " & " + str(mut_freq) + "\n")
                print(protocol, " & ", round(statistics.mean(results[protocol]["ANI"]), 2), " & ", round(statistics.median(results[protocol]["ANI"]), 2), " & ", mut_freq)
            else:
                for params in results[protocol]:
                    log_file.write(protocol + " & " + str(params) + " & " + ", ".join([str(round(r, 4)) for r in results[protocol][params]["ANI"]]) + " & " + str(mut_freq) + "\n")
                    print(protocol, " & ", params, " & ", round(statistics.mean(results[protocol][params]["ANI"]), 2), " & ", round(statistics.median(results[protocol][params]["ANI"]),2), " & ", mut_freq)

    log_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Estimate ANI on simulated sequences using Strobemer-adapted Mash distance", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--L', type=int, default=1000, help='Length of simulated sequences')
    parser.add_argument('--nr_exp', type=int, default=1000, help='Number of simulated experiments')
    parser.add_argument('--experiment_type', type=str, default="all", help='experiment type choose between "all", "controlled", "specified" or "only_subs"')
    parser.add_argument('--subs_freq', type=float, default=0.33, help='substitution frequency among all mutations for --experiment_type "specified"; rest split evenly in insertions and deletions')
    parser.add_argument('--mut_freqs', nargs='+', type=float, default=[0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10], help='mutation frequencies [0,1]')
    parser.add_argument('--k_size', type=int, default=30, help='k-mer/strobemer length')
    parser.add_argument('--w', type=int, default=1, help='number of hashes used in a sliding window for thinning (w=1 means no thinning)')
    #parser.add_argument('--orders', type=list, default=[4, ], help='List with orders of strobes to be analzyed')
    parser.add_argument('--orders', nargs='+', type=int, default=[2, ], help='List with orders of strobes to be analzyed')
    parser.add_argument('--w_low', type=int, default=25, help='minimum window offset to the previous window (wMin > 0)')
    parser.add_argument('--w_high', type=int, default=50, help='maximum window offset to the previous window (wMin <= wMax)')
    parser.add_argument('--strobe_fractions', nargs='+', type=float, default=[0.8], help='Fraction of sampled strobemers, rest kmers')
    parser.add_argument('--all_methods', action="store_true", help='perform matching analysis on simulated data for all (mixed-)strobemer and (mixed-)altstrobe seeding techniques')
    parser.add_argument('--method', type=str, default="none", help='choose seeding technique')
    parser.add_argument('--altstrobes_generalized', action="store_true", help='perform matching analysis on simulated data for altstrobes of all combinations from (1,k-1) to (k/2,k/2)')
    parser.add_argument('--altstrobes_size_distribution', action="store_true", help='perform matching analysis on simulated data for altstrobes with strobe size distribution (k_s, k_l) determined by strobe_fraction')
    parser.add_argument('--multistrobes', action="store_true", help='perform matching analysis on simulated data for multistrobes of all combinations from (k_boundary,k-k_boundary) to (k/2,k/2)')
    parser.add_argument('--multistrobes_size_distribution', action="store_true", help='perform matching analysis on simulated data for multistrobes of all combinations from (k_boundary,k-k_boundary) to (k/2,k/2) with strobe size distribution (k_s, k_l) determined by strobe_fraction')
    parser.add_argument('--k_boundary', type=int, default=5, help='minimum strobe length (k >= 4 recommended to ensure uniqueness)')
    parser.add_argument('--mixedstrobes', action="store_true", help='perform matching analysis on simulated data for user defined mixed seeding techniques')
    parser.add_argument('--mixedstrobes_methods', type=list, default=["randstrobes", "kmers"], help='List with two seeding methods to sample mixedstrobes')
    parser.add_argument('--verbose', action="store_true")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    # methods = (
    #     "kmers", "spaced_kmers_dense", "spaced_kmers_sparse",
    #     "minstrobes", "randstrobes", "hybridstrobes", "altstrobes",
    #     "mixedminstrobes", "mixedrandstrobes", "mixedhybridstrobes", "mixedaltstrobes")
    methods = (
        "kmers", "spaced_kmers_dense", "spaced_kmers_sparse",
        "minstrobes", "randstrobes", "hybridstrobes", "altstrobes",
        "mixedrandstrobes", "multistrobes")

    if args.altstrobes_generalized:
        args.methods = ("altstrobes_generalized",)
        args.strobe_lengths = [(k1, args.k_size-k1) for k1 in range(1, int(args.k_size/2)+1)]
    elif args.multistrobes:
        args.methods = ("multistrobes",)
    elif args.altstrobes_size_distribution:
        args.methods = ("altstrobes_size_distribution",)
    elif args.multistrobes_size_distribution:
        args.methods = ("multistrobes_size_distribution",)
    elif args.mixedstrobes:
        args.methods = ("mixedstrobes",)
        args.method1 = args.mixedstrobes_methods[0]
        args.method2 = args.mixedstrobes_methods[1]
    elif args.method != "none":
        assert args.method in methods, "[Error] Seeding technique not implemented"
        args.methods = (args.method,)
    else:
        args.methods = methods

    main(args)
