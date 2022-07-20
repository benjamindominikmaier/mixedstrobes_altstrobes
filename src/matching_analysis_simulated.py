#!/usr/bin/env python3.9
# -*- coding: utf-8 -*
import os
import sys
import argparse
import random
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot

from modules import help_functions
from modules import indexing_Maier_altstrobes as indexing


def get_match_coverage(seq_len: int, mers: dict, matches: dict, order: int,
                       span: int) -> int:
    """
    Header

    :param seq_len:
    :param mers:
    :param matches:
    :param order: number of substrings/strobes
    :param span: Span is k if kmer, length between first and last position if spaced kmer, and span between the start of the first strobe and the last nucleotide of the last strobe.
    :returns:
    """
    if not matches:
        return 0

    if order == 1:
        match_list = sorted([(p, p+span) for (p, h) in mers.items() if h in matches])
    else:
        match_list = sorted([(p[0], p[-1]+span) for (p, h) in mers.items() if h in matches])

    covered_bases = match_list[0][1] - match_list[0][0]
    max_stop = match_list[0][1]
    for start, stop in match_list:
        if start < max_stop and stop > max_stop:
            covered_bases += stop - max_stop
            max_stop = stop
        elif start >= max_stop:
            covered_bases += stop - start
            max_stop = stop
    # print(covered_bases)
    return covered_bases  # round(100*covered_bases/seq_len, 1)


def seq_covered_spaced_kmers(mers: dict, matches: dict, seq: str, positions) -> int:
    """
    Function specific to calculate the coverage of spaced k-mers
    since we have all the "sampled" positions in a spaced k-mer
    we can keep an active set of covered positions as we iterate
    through the string.

    :param mers:
    :param matches:
    :param seq: a string with a nucleotide sequence
    :param positions:
    :returns:
    """
    seq_covered = 0
    if not matches:
        return seq_covered

    spaced_kmer_sampled_pos = sorted(positions)
    covered_pos = [0]*len(seq)
    all_match_pos_vector = sorted([p for (p, k) in mers.items() if k in matches])

    for p in all_match_pos_vector:
        for j in spaced_kmer_sampled_pos:
            covered_pos[p + j] += 1

    # all positions covered at least once
    c = sum([1 for p in covered_pos if p > 0])
    return c


def get_sequence_coverage(mers: dict, matches: dict, order: int, k_len: int) -> int:
    """
    Header

    :param mers:
    :param matches:
    :param order: number of substrings/strobes
    :param k_len:
    :returns:
    """
    covered_bases = 0
    if not matches:
        return covered_bases

    if order == 1:
        all_pos_vector = sorted([p for (p, k) in mers.items() if k in matches])
    else:
        match_list = [list(p) for (p, k) in mers.items() if k in matches]
        all_pos_vector = sorted([p for sublist in match_list for p in sublist])

    prev_p = all_pos_vector[0]
    covered_bases += k_len  # for first pos
    if len(all_pos_vector) == 1:
        return covered_bases

    for p in all_pos_vector[1:]:
        if p <= prev_p + k_len - 1:
            covered_bases += p-prev_p
        else:
            covered_bases += k_len

        prev_p = p

    return covered_bases


def get_intervals(mers: dict, matches: dict, order: int) -> tuple:
    """
    Header

    :param mers:
    :param matches:
    :param order: number of substrings/strobes
    :returns:
    """
    if not matches:
        return [], []

    if order == 1:
        all_pos_vector = sorted([p for (p, k) in mers.items() if k in matches])
    else:
        match_list = [[p for p in range(pos[0], pos[-1]+1)] for (pos, k) in mers.items() if k in matches]
        all_pos_vector = sorted(set([p for sublist in match_list for p in sublist]))

    ivls = []
    iv_start = all_pos_vector[0]
    length = 0
    for p1, p2 in zip(all_pos_vector[:-1], all_pos_vector[1:]):
        if p2 == p1 + 1:
            length += 1
        # elif p2 == p1: # for the strobes
        #     pass
        elif p2 > p1 + 1:
            ivls.append((iv_start, iv_start+length))
            length = 0
            iv_start = p2

    if len(all_pos_vector) > 1:
        if p2 <= p1 + 1:
            ivls.append((iv_start, iv_start+length))
    elif len(all_pos_vector) == 1:
        ivls.append((iv_start, iv_start))
    # print(ivls)
    return ivls, all_pos_vector


def statistics(ivls, seq: str, k: int) -> tuple:
    """
    Header

    :param ivls:
    :param seq: a string with a nucleotide sequence
    :param k:
    :returns:
    """
    if not ivls:
        return 1, [len(seq)], 0
    seq_covered = 0
    nr_islands = 0
    gap_lengths = []

    prev_stop = 0  # ivls[0][1] + k
    for i, (start, stop) in enumerate(ivls):
        if i == 0:
            seq_covered += (stop-start) + k
            if start > 0:
                nr_islands += 1
                gap_lengths.append(start - prev_stop)

        elif start > prev_stop + k:
            seq_covered += (stop-start) + k
            nr_islands += 1
            # Extra check to see if the gap is at least 1 nt. We may end up here because
            # it may be that start = prev_stop + k + 1, which leads to a 0nt gap
            # (0nt can occur because of homopolymer stretches)
            # if start > prev_stop + k + 1:
            gap_lengths.append(start - prev_stop - k)
            assert start - prev_stop >= 0, "found: {0},{1},{2}: {3}".format(start, prev_stop, k, ivls)
        else:  # overlapping previous hit
            seq_covered += (stop-start) + k - (k - (start - prev_stop))

        prev_stop = stop

    if ivls[-1][0] + k - 1 < len(seq):
        nr_islands += 1
        gap_lengths.append(len(seq) - prev_stop - k + 1)

    return nr_islands, gap_lengths, seq_covered


def analyze_strobemers(seq1: str, seq2: str, k_size: int, order: int,
                       hash_fcn: str, w: int, w_low: int = 0, w_high: int = 50,
                       fraction: float = 0.5) -> tuple:
    """
    Header

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns:
    """
    assert k_size % order == 0, "Not even kmer length, results will be different"
    if hash_fcn in ("mixedminstrobes", "mixedminstrobes_generalized", "mixedrandstrobes", "mixedrandstrobes_generalized", "mixedhybridstrobes", "mixedhybridstrobes_generalized"):
        strobemers1 = getattr(indexing, hash_fcn)(seq1, k_size, w_low, w_high, w, order, fraction)
        strobemers2 = getattr(indexing, hash_fcn)(seq2, k_size, w_low, w_high, w, order, fraction)
    else:
        strobemers1 = getattr(indexing, hash_fcn)(seq1, k_size, w_low, w_high, w, order)
        strobemers2 = getattr(indexing, hash_fcn)(seq2, k_size, w_low, w_high, w, order)

    matches = set(strobemers1.values()) & set(strobemers2.values())
    m = len(matches)
    mp = len(strobemers1.values())
    ivls, all_pos_vector = get_intervals(strobemers1, matches, order)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size//order)
    seq_cov = get_sequence_coverage(strobemers1, matches, order, k_size//order)
    match_coverage = get_match_coverage(len(seq1), strobemers1, matches, order, k_size//order)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


def analyze_kmers(seq1: str, seq2: str, k_size: int, w: int) -> tuple:
    """
    Header

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of the kmers
    :param w: number of kmers used in a sliding window for thinning (w=1 means no thinning)
    :returns:
    """
    # kmers
    kmers_pos1 = indexing.kmers(seq1, k_size, w)
    kmers_pos2 = indexing.kmers(seq2, k_size, w)

    matches = set(kmers_pos1.values()) & set(kmers_pos2.values())
    m = len(matches)
    mp = len(kmers_pos1.values())
    ivls, all_pos_vector = get_intervals(kmers_pos1, matches, 1)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size)
    seq_cov = get_sequence_coverage(kmers_pos1, matches, 1, k_size)
    match_coverage = get_match_coverage(len(seq1), kmers_pos1, matches, 1, k_size)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


def analyze_spaced_kmers(seq1: str, seq2: str, k_size: int, span_size: int, w: int) -> tuple:
    """
    Header

    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of the kmers
    :param span_size: length between first and last position
    :param w: number of spaced_kmers used in a sliding window for thinning (w=1 means no thinning)
    :returns:
    """
    positions = set(random.sample(range(1, span_size-1), k_size-2))
    positions.add(0)
    positions.add(span_size - 1)  # asserts first and last position is sampled so that we have a spaced kmer of length span size
    spaced_kmers_seq1 = indexing.spaced_kmers(seq1, k_size, span_size, positions, w)
    spaced_kmers_seq2 = indexing.spaced_kmers(seq2, k_size, span_size, positions, w)
    matches = set(spaced_kmers_seq1.values()) & set(spaced_kmers_seq2.values())
    m = len(matches)
    mp = len(spaced_kmers_seq1.values())
    ivls, all_pos_vector = get_intervals(spaced_kmers_seq1, matches, 1)
    nr_islands, gap_lengths, _ = statistics(ivls, seq1, span_size)
    # we compute coverage for spaced k-mers with specific function
    seq_cov = seq_covered_spaced_kmers(spaced_kmers_seq1, matches, seq1, positions)
    match_coverage = get_match_coverage(len(seq1), spaced_kmers_seq1, matches, 1, span_size)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


def analyze_altstrobes(seq1: str, seq2: str, k_size: int, order: int,
                       hash_fcn: str, w: int, w_low: int = 0, w_high: int = 50,
                       fraction: float = 0.5) -> tuple:
    """
    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns:
    """
    assert k_size % order == 0, "Not even kmer length, results will be different"
    altstrobes1 = indexing.altstrobes(seq1, k_size, w_low, w_high, w, order)
    altstrobes2 = indexing.altstrobes(seq2, k_size, w_low, w_high, w, order)

    order = order+1

    matches = set(altstrobes1.values()) & set(altstrobes2.values())
    m = len(matches)
    mp = len(altstrobes1.values())
    ivls, all_pos_vector = get_intervals(altstrobes1, matches, order)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size//order)
    seq_cov = get_sequence_coverage(altstrobes1, matches, order, k_size//order+1)
    match_coverage = get_match_coverage(len(seq1), altstrobes1, matches, order, k_size//order)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


def analyze_mixedaltstrobes(seq1: str, seq2: str, k_size: int, order: int,
                       hash_fcn: str, w: int, w_low: int = 0, w_high: int = 50,
                       fraction: float = 0.5) -> tuple:
    """
    :param seq1: a string with a nucleotide sequence which is compared to seq2
    :param seq2: a string with a nucleotide sequence which is compared to seq1
    :param k_size: length of all strobes/substrings (len(strobe_1) +  ... + len(strobe_n))
    :param order: number of substrings/strobes
    :param hash_fcn: a string with the function name of the kmer/strobemer protocol
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns:
    """
    assert k_size % order == 0, "Not even kmer length, results will be different"
    altstrobes1 = indexing.mixedaltstrobes(seq1, k_size, w_low, w_high, w, order, fraction)
    altstrobes2 = indexing.mixedaltstrobes(seq2, k_size, w_low, w_high, w, order, fraction)

    order = order+1

    matches = set(altstrobes1.values()) & set(altstrobes2.values())
    m = len(matches)
    mp = len(altstrobes1.values())
    ivls, all_pos_vector = get_intervals(altstrobes1, matches, order)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size//order)
    seq_cov = get_sequence_coverage(altstrobes1, matches, order, k_size//order)
    match_coverage = get_match_coverage(len(seq1), altstrobes1, matches, order, k_size//order)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


def print_matches(all_pos_vector, method: str):
    """
    Header

    :param all_pos_vector:
    :param method:
    :returns:
    """
    s = set(all_pos_vector)
    for i in range(100):
        if i in s:
            print("X", end='')
        else:
            print(" ", end='')
    print(method)


def get_e_size(all_islands, L: int, nr_exp: int) -> float:
    """
    Header

    :param all_islands:
    :param L: an integer representing the desired sequence length for the analysis
    :param nr_exp: in integer representing the number of desired experiments
    :returns:
    """
    # print("all_islands",all_islands)
    sum_of_squares = sum([x**2 for x in all_islands])
    return sum_of_squares/(L*nr_exp)


def main(args):
    """
    """
    L = 10000 # 10000
    k_size = 30
    nr_exp = 100  # 1000
    w = 1  # thinning, w = 1  means no thinning. w =1, 10, 20 was used in the evaluations.
    mut_freqs = [0.01, 0.05, 0.1]  # [0.1]
    orders = [2, ]
    w_low = 25
    w_high = 50
    methods = ("kmers", "spaced_kmers_dense", "spaced_kmers_sparse", "minstrobes", "randstrobes", "hybridstrobes", "altstrobes", "mixedminstrobes", "mixedrandstrobes", "mixedhybridstrobes", "mixedaltstrobes")
    mixedstrobe_fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # experiment_type choose between 'only_subs', 'controlled' or 'all'
    experiment_type = "all"  # "controlled" # "all" #"only_subs" # "" # for spaced kmers
    # mut_freq = 0.5 #0.01 #, 0.05, 0.1]
    list_for_illustration = [[], [], [], [], [], []]

    for mut_freq in mut_freqs:
        print("MUTATION RATE:", mut_freq)
        results = dict()

        for exp_id in range(nr_exp):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])

            # controlled or random experiment
            if experiment_type == 'only_subs':
                muts = set(random.sample(range(len(seq1)), int(L*mut_freq)))
                seq2 = "".join([
                    seq1[i] if i not in muts
                    else random.choice([help_functions.reverse_complement(seq1[i])])
                    for i in range(len(seq1))
                ])
            elif experiment_type == 'controlled':
                # muts = set(range(15,L,15)) # every 15th nt for figure 2 only!
                muts = set(range(20, L, 20))
                seq2 = "".join([
                    seq1[i] if i not in muts
                    else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")])
                    for i in range(len(seq1))
                ])
            elif experiment_type == 'all':
                muts = set(random.sample(range(len(seq1)), int(L*mut_freq)))
                seq2 = "".join([
                    seq1[i] if i not in muts
                    else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")])
                    for i in range(len(seq1))
                ])
            else:
                print("Wrong experiment label specified")
                sys.exit()

            for hash_fcn in methods:
                if hash_fcn == "kmers":
                    results.setdefault("kmers", {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                    m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_kmers(seq1, seq2, k_size, w)
                    results["kmers"]["m"] += m
                    results["kmers"]["sc"] += sc
                    results["kmers"]["gaps"].append(gaps)
                    results["kmers"]["mc"] += match_coverage
                    results["kmers"]["mp"] += mp

                elif hash_fcn == "spaced_kmers_dense":
                    results.setdefault("spaced_kmers_dense", {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                    m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, k_size, k_size+k_size//2, w)
                    results["spaced_kmers_dense"]["m"] += m
                    results["spaced_kmers_dense"]["sc"] += sc
                    results["spaced_kmers_dense"]["gaps"].append(gaps)
                    results["spaced_kmers_dense"]["mc"] += match_coverage
                    results["spaced_kmers_dense"]["mp"] += mp

                elif hash_fcn == "spaced_kmers_sparse":
                    results.setdefault("spaced_kmers_sparse", {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                    m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, k_size, 3*k_size, w)
                    results["spaced_kmers_sparse"]["m"] += m
                    results["spaced_kmers_sparse"]["sc"] += sc
                    results["spaced_kmers_sparse"]["gaps"].append(gaps)
                    results["spaced_kmers_sparse"]["mc"] += match_coverage
                    results["spaced_kmers_sparse"]["mp"] += mp

                elif hash_fcn == "altstrobes":
                        results.setdefault("altstrobes", {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                        m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_altstrobes(seq1, seq2, k_size, 2, hash_fcn, w, w_low=w_low, w_high=w_high)
                        results["altstrobes"]["m"] += m
                        results["altstrobes"]["sc"] += sc
                        results["altstrobes"]["gaps"].append(gaps)
                        results["altstrobes"]["mc"] += match_coverage
                        results["altstrobes"]["mp"] += mp

                elif hash_fcn == "mixedaltstrobes":
                    results.setdefault(hash_fcn, dict())
                    for order in orders:
                        for fraction in mixedstrobe_fractions:
                            results["mixedaltstrobes"].setdefault((order, k_size//order, w_low, w_high, fraction), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_mixedaltstrobes(seq1, seq2, k_size, order, hash_fcn, w, w_low=w_low, w_high=w_high, fraction=fraction)
                            results["mixedaltstrobes"][(order, k_size//order, w_low, w_high, fraction)]["m"] += m
                            results["mixedaltstrobes"][(order, k_size//order, w_low, w_high, fraction)]["sc"] += sc
                            results["mixedaltstrobes"][(order, k_size//order, w_low, w_high, fraction)]["gaps"].append(gaps)
                            results["mixedaltstrobes"][(order, k_size//order, w_low, w_high, fraction)]["mc"] += match_coverage
                            results["mixedaltstrobes"][(order, k_size//order, w_low, w_high, fraction)]["mp"] += mp

                elif hash_fcn in ("mixedminstrobes", "mixedminstrobes_generalized", "mixedrandstrobes", "mixedrandstrobes_generalized", "mixedhybridstrobes", "mixedhybridstrobes_generalized"):
                    results.setdefault(hash_fcn, dict())
                    for order in orders:
                        for fraction in mixedstrobe_fractions:
                            results[hash_fcn].setdefault((order, k_size//order, w_low, w_high, fraction), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, order, hash_fcn, w, w_low=w_low, w_high=w_high, fraction=fraction)
                            results[hash_fcn][(order, k_size//order, w_low, w_high, fraction)]["m"] += m
                            results[hash_fcn][(order, k_size//order, w_low, w_high, fraction)]["sc"] += sc
                            results[hash_fcn][(order, k_size//order, w_low, w_high, fraction)]["gaps"].append(gaps)
                            results[hash_fcn][(order, k_size//order, w_low, w_high, fraction)]["mc"] += match_coverage
                            results[hash_fcn][(order, k_size//order, w_low, w_high, fraction)]["mp"] += mp

                else:
                    results.setdefault(hash_fcn, dict())
                    for order in orders:
                        results[hash_fcn].setdefault((order, k_size//order, w_low, w_high), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                        m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, order, hash_fcn, w, w_low=w_low, w_high=w_high)
                        results[hash_fcn][(order, k_size//order, w_low, w_high)]["m"] += m
                        results[hash_fcn][(order, k_size//order, w_low, w_high)]["sc"] += sc
                        results[hash_fcn][(order, k_size//order, w_low, w_high)]["gaps"].append(gaps)
                        results[hash_fcn][(order, k_size//order, w_low, w_high)]["mc"] += match_coverage
                        results[hash_fcn][(order, k_size//order, w_low, w_high)]["mp"] += mp

        for protocol in results:
            if protocol == "kmers" or protocol == "spaced_kmers_sparse" or protocol == "spaced_kmers_dense" or protocol == "altstrobes":
                flat = [g for l in results[protocol]["gaps"] for g in l]
                if flat:
                    # avg_island_len = sum(flat)/len(flat)
                    # print(protocol)
                    e_size = get_e_size(flat, L, nr_exp)
                # else:
                #     avg_island_len = 0
                res = [
                    round(100*results[protocol]["m"]/results[protocol]["mp"], 1),
                    100*results[protocol]["sc"]/(L*nr_exp),
                    100*results[protocol]["mc"]/(L*nr_exp),
                    e_size
                ]
                print(protocol, " & ".join([str(round(r, 1)) for r in res]))
            else:
                for params in results[protocol]:
                    flat = [g for l in results[protocol][params]["gaps"] for g in l]
                    if flat:
                        # avg_island_len = sum(flat)/len(flat)
                        # print(protocol, params)
                        e_size = get_e_size(flat, L, nr_exp)
                    # else:
                        # avg_island_len = 0
                    res = [
                        round(100*results[protocol][params]["m"]/results[protocol][params]["mp"], 1),
                        100*results[protocol][params]["sc"]/(L*nr_exp),
                        100*results[protocol][params]["mc"]/(L*nr_exp),
                        e_size
                    ]
                    print(protocol, params, " & ".join([str(round(r, 1)) for r in res]))

    # print(results)

    # # random mutation mositions
    # for mut_freq in [0.01, 0.05, 0.1]:
    #     for exp_id in range(10):
    #         seq1 = "".join([random.choice("ACGT") for i in range(L)])
    #         muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
    #         # muts = set(range(20,1000,20)) #set([20,40,60,80])
    #         # print(muts)
    #         seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
    #         print()
    #         print("MUT FREQ:", mut_freq)
    #         positions_matching_kmers(seq1, seq2, k)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    # parser.add_argument('--k', type=int, default=13, help='Kmer size')
    # parser.add_argument('--w', type=int, default=20, help='Window size')
    parser.add_argument('--outfolder', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()

    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit()

    main(args)
