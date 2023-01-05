#!/usr/bin/env python3.9
# -*- coding: utf-8 -*
import os
import sys
import argparse
import random
import numpy as np
import collections
from typing import Iterator

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot
from itertools import chain

from modules import help_functions
from modules import indexing_Maier_altstrobes as indexing


def get_match_coverage(seq_len: int, mers: dict, matches: dict, order: int,
                       span: int) -> int:
    """
    Computes the proportion of nucleotides covered by the k-mers and strobemers
    from end-to-end including potential gaps

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
    Computes the proportion of nucleotides covered by the strobes of matches,
    this function distinguishes from match coverage by disregarding the gaps between the strobes

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


def gap_finder(gap_pos):
    """
    """
    count = 0
    gap_lengths = []
    for i in range(len(gap_pos) - 1):
        # Check if the next number is consecutive
        if gap_pos[i] + 1 == gap_pos[i+1]:
            count += 1
        else:
            # If it is not append the count and restart counting
            gap_lengths.append(count)
            count = 1
    # Since we stopped the loop one early append the last count
    gap_lengths.append(count)
    return gap_lengths


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

    matches = set(strobemers1.values()) & set(strobemers2.values())
    m = len(matches)
    mp = len(strobemers1.values())
    ivls, all_pos_vector = get_intervals(strobemers1, matches, order)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size//order)
    seq_cov = get_sequence_coverage(strobemers1, matches, order, k_size//order)
    match_coverage = get_match_coverage(len(seq1), strobemers1, matches, order, k_size//order)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


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

    matches = set(altstrobes1.values()) & (set(altstrobes2.values()) | set(altstrobes2_rev.values()))

    m = len(matches)
    mp = len(altstrobes1.values())

    hash_values_tmp = list(altstrobes1.values())
    hash_values = collections.Counter(hash_values_tmp)

    altstrobe_positions = [key for key, value in altstrobes1.items() if value in matches]

    sc_pos = []
    mc_pos = []

    for x1, x2, x3 in altstrobe_positions:
        if x1 == x2:  # first k2, then k1
            sc_pos.append(range(x1, x1+k_size2))
            sc_pos.append(range(x3, x3+k_size1))
            mc_pos.append(range(x1, x3+k_size1))
        else:
            sc_pos.append(range(x1, x1+k_size1))
            sc_pos.append(range(x3, x3+k_size2))
            mc_pos.append(range(x1, x3+k_size2))

    sc_pos = set(chain.from_iterable(sc_pos))
    mc_pos = set(chain.from_iterable(mc_pos))
    gap_pos = [pos for pos in range(len(seq1)) if pos not in mc_pos]
    gaps = gap_finder(gap_pos)

    return m, mp, len(sc_pos), len(mc_pos), gaps


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

    matches = set(multistrobes1.values()) & (set(multistrobes2.values()))

    m = len(matches)
    mp = len(multistrobes1.values())

    hash_values_tmp = list(multistrobes1.values())
    hash_values = collections.Counter(hash_values_tmp)

    multistrobe_positions = [key for key, value in multistrobes1.items() if value in matches]

    # stat = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0}
    # for x in multistrobe_positions:
    #     s1 = len(x[0])
    #     s2 = len(x[1])
    #     assert order*(s1 + s2)/2 == k_size; "[Error] " + str(s1) + str(s2)
    #     stat[s1] += 1
    # print(stat)

    sc_pos = []
    mc_pos = []

    for x in multistrobe_positions:
        for elem in x:
            sc_pos.append(elem)
        mc_pos.append(range(min(x[0]), max(x[-1])+1))


    sc_pos = set(chain.from_iterable(sc_pos))
    mc_pos = set(chain.from_iterable(mc_pos))
    gap_pos = [pos for pos in range(len(seq1)) if pos not in mc_pos]
    gaps = gap_finder(gap_pos)

    return m, mp, len(sc_pos), len(mc_pos), gaps


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

    matches = set(altstrobes1.values()) & set(altstrobes2.values())
    m = len(matches)
    mp = len(altstrobes1.values())
    ivls, all_pos_vector = get_intervals(altstrobes1, matches, order)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size//order)
    seq_cov = get_sequence_coverage(altstrobes1, matches, order, k_size//order)
    match_coverage = get_match_coverage(len(seq1), altstrobes1, matches, order, k_size//order)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


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

    matches = set(altstrobes1.values()) & set(altstrobes2.values())
    m = len(matches)
    mp = len(altstrobes1.values())
    ivls, all_pos_vector = get_intervals(altstrobes1, matches, order)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size//order)
    seq_cov = get_sequence_coverage(altstrobes1, matches, order, k_size//order)
    match_coverage = get_match_coverage(len(seq1), altstrobes1, matches, order, k_size//order)

    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


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

    matches = set(mixedstrobes1.values()) & set(mixedstrobes2.values())
    m = len(matches)
    mp = len(mixedstrobes1.values())
    mixedstrobe_positions = [key for key, value in mixedstrobes1.items() if value in matches]

    sc_pos = []
    mc_pos = []

    for x1, x2 in mixedstrobe_positions:
        sc_pos.append(x1)
        sc_pos.append(x2)
        mc_pos.append(range(min(x1), max(x2)+1))


    sc_pos = set(chain.from_iterable(sc_pos))
    mc_pos = set(chain.from_iterable(mc_pos))
    gap_pos = [pos for pos in range(len(seq1)) if pos not in mc_pos]
    gaps = gap_finder(gap_pos)

    return m, mp, len(sc_pos), len(mc_pos), gaps


def get_e_size(all_islands: list, L: int, nr_exp: int) -> float:
    """
    Computes the expected island size, whereby an island is the maximal interval
    of consecutive nucleotides without matches:

    :param all_islands: list with island size lengths
    :param L: an integer representing the desired sequence length for the analysis
    :param nr_exp: in integer representing the number of desired experiments
    :returns: float indicating the expected island size
    """
    # print("all_islands",all_islands)
    sum_of_squares = sum([x**2 for x in all_islands])
    return sum_of_squares/(L*nr_exp)


def main(args):
    """
    Computes matching metrics (m, mp, sc, mc, e-size) for simulated sequences
    using various seeding techniques (kmer-based, strobemers, altstrobes)
    """

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
                    results.setdefault("kmers", {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                    m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_kmers(seq1, seq2, args.k_size, args.w)
                    results["kmers"]["m"] += m
                    results["kmers"]["sc"] += sc
                    results["kmers"]["gaps"].append(gaps)
                    results["kmers"]["mc"] += match_coverage
                    results["kmers"]["mp"] += mp

                elif hash_fcn == "spaced_kmers_dense":
                    results.setdefault("spaced_kmers_dense", {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                    m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, args.k_size, args.k_size+args.k_size//2, args.w)
                    results["spaced_kmers_dense"]["m"] += m
                    results["spaced_kmers_dense"]["sc"] += sc
                    results["spaced_kmers_dense"]["gaps"].append(gaps)
                    results["spaced_kmers_dense"]["mc"] += match_coverage
                    results["spaced_kmers_dense"]["mp"] += mp

                elif hash_fcn == "spaced_kmers_sparse":
                    results.setdefault("spaced_kmers_sparse", {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                    m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, args.k_size, 3*args.k_size, args.w)
                    results["spaced_kmers_sparse"]["m"] += m
                    results["spaced_kmers_sparse"]["sc"] += sc
                    results["spaced_kmers_sparse"]["gaps"].append(gaps)
                    results["spaced_kmers_sparse"]["mc"] += match_coverage
                    results["spaced_kmers_sparse"]["mp"] += mp

                elif hash_fcn == "altstrobes_size_distribution":
                    results.setdefault("altstrobes_size_distribution", dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results["altstrobes_size_distribution"].setdefault((order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_altstrobes(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, arg=100*fraction)
                            results["altstrobes_size_distribution"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["m"] += m
                            results["altstrobes_size_distribution"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["sc"] += sc
                            results["altstrobes_size_distribution"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["gaps"].append(gaps)
                            results["altstrobes_size_distribution"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["mc"] += match_coverage
                            results["altstrobes_size_distribution"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["mp"] += mp

                elif hash_fcn == "altstrobes":
                    results.setdefault("altstrobes", dict())
                    for order in args.orders:
                        results["altstrobes"].setdefault((order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                        m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_altstrobes(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high)
                        results["altstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high)]["m"] += m
                        results["altstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high)]["sc"] += sc
                        results["altstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high)]["gaps"].append(gaps)
                        results["altstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high)]["mc"] += match_coverage
                        results["altstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high)]["mp"] += mp

                elif hash_fcn == "altstrobes_generalized":
                    results.setdefault("altstrobes_generalized", dict())
                    for k1, k2 in args.strobe_lengths:
                        results["altstrobes_generalized"].setdefault((k1, k2), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                        m, mp, sc, match_coverage, gaps = analyze_altstrobes_generalized(seq1, seq2, k1, k2, 2, args.w, w_low=args.w_low, w_high=args.w_high)
                        results["altstrobes_generalized"][(k1, k2)]["m"] += m
                        results["altstrobes_generalized"][(k1, k2)]["sc"] += sc
                        results["altstrobes_generalized"][(k1, k2)]["gaps"].append(gaps)
                        results["altstrobes_generalized"][(k1, k2)]["mc"] += match_coverage
                        results["altstrobes_generalized"][(k1, k2)]["mp"] += mp

                elif hash_fcn == "mixedaltstrobes":
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results["mixedaltstrobes"].setdefault((order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_mixedaltstrobes(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, fraction=fraction)
                            results["mixedaltstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["m"] += m
                            results["mixedaltstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["sc"] += sc
                            results["mixedaltstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["gaps"].append(gaps)
                            results["mixedaltstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["mc"] += match_coverage
                            results["mixedaltstrobes"][(order, (2*args.k_size)//(3*order), (4*args.k_size)//(3*order), args.w_low, args.w_high, fraction)]["mp"] += mp

                elif hash_fcn == "multistrobes":
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        results["multistrobes"].setdefault((order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                        m, mp, sc, match_coverage, gaps = analyze_multistrobes(seq1, seq2, args.k_size, order, args.w, w_low=args.w_low, w_high=args.w_high, k_boundary=args.k_boundary)
                        results["multistrobes"][(order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1)]["m"] += m
                        results["multistrobes"][(order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1)]["sc"] += sc
                        results["multistrobes"][(order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1)]["gaps"].append(gaps)
                        results["multistrobes"][(order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1)]["mc"] += match_coverage
                        results["multistrobes"][(order, args.k_size, args.k_boundary, int(2*args.k_size/order-args.k_boundary), args.w_low, args.w_high, 1)]["mp"] += mp

                elif hash_fcn == "multistrobes_size_distribution":
                    results.setdefault("multistrobes_size_distribution", dict())
                    for fraction in args.strobe_fractions:
                        results["multistrobes_size_distribution"].setdefault((2, args.k_size, args.w_low, args.w_high, fraction), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                        m, mp, sc, match_coverage, gaps = analyze_multistrobes(seq1, seq2, args.k_size, 2, args.w, w_low=args.w_low, w_high=args.w_high, k_boundary=args.k_boundary, arg=100*fraction)
                        results["multistrobes_size_distribution"][(2, args.k_size, args.w_low, args.w_high, fraction)]["m"] += m
                        results["multistrobes_size_distribution"][(2, args.k_size, args.w_low, args.w_high, fraction)]["sc"] += sc
                        results["multistrobes_size_distribution"][(2, args.k_size, args.w_low, args.w_high, fraction)]["gaps"].append(gaps)
                        results["multistrobes_size_distribution"][(2, args.k_size, args.w_low, args.w_high, fraction)]["mc"] += match_coverage
                        results["multistrobes_size_distribution"][(2, args.k_size, args.w_low, args.w_high, fraction)]["mp"] += mp

                elif hash_fcn == "mixedstrobes":
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results["mixedstrobes"].setdefault((args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                            m, mp, sc, match_coverage, gaps = analyze_mixedstrobes(args.method1, args.method2, seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, fraction=fraction)
                            results["mixedstrobes"][(args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction)]["m"] += m
                            results["mixedstrobes"][(args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction)]["sc"] += sc
                            results["mixedstrobes"][(args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction)]["gaps"].append(gaps)
                            results["mixedstrobes"][(args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction)]["mc"] += match_coverage
                            results["mixedstrobes"][(args.method1, args.method2, order, args.k_size//order, args.w_low, args.w_high, fraction)]["mp"] += mp

                elif hash_fcn in ("mixedminstrobes", "mixedminstrobes_generalized", "mixedrandstrobes", "mixedrandstrobes_generalized", "mixedhybridstrobes", "mixedhybridstrobes_generalized"):
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        for fraction in args.strobe_fractions:
                            results[hash_fcn].setdefault((order, args.k_size//order, args.w_low, args.w_high, fraction), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high, fraction=fraction)
                            results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high, fraction)]["m"] += m
                            results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high, fraction)]["sc"] += sc
                            results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high, fraction)]["gaps"].append(gaps)
                            results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high, fraction)]["mc"] += match_coverage
                            results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high, fraction)]["mp"] += mp

                else:
                    results.setdefault(hash_fcn, dict())
                    for order in args.orders:
                        results[hash_fcn].setdefault((order, args.k_size//order, args.w_low, args.w_high), {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0})
                        m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, args.k_size, order, hash_fcn, args.w, w_low=args.w_low, w_high=args.w_high)
                        results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high)]["m"] += m
                        results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high)]["sc"] += sc
                        results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high)]["gaps"].append(gaps)
                        results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high)]["mc"] += match_coverage
                        results[hash_fcn][(order, args.k_size//order, args.w_low, args.w_high)]["mp"] += mp

        for protocol in results:
            if protocol == "kmers" or protocol == "spaced_kmers_sparse" or protocol == "spaced_kmers_dense":
                flat = [g for l in results[protocol]["gaps"] for g in l]
                if flat:
                    # avg_island_len = sum(flat)/len(flat)
                    # print(protocol)
                    e_size = get_e_size(flat, args.L, args.nr_exp)
                # else:
                #     avg_island_len = 0
                res = [
                    round(100*results[protocol]["m"]/results[protocol]["mp"], 1),
                    100*results[protocol]["sc"]/(args.L*args.nr_exp),
                    100*results[protocol]["mc"]/(args.L*args.nr_exp),
                    e_size
                ]
                print(protocol, " & ", args.orders[0], " & - & 0 &", " & ".join([str(round(r, 1)) for r in res]), " & ", mut_freq)
            else:
                for params in results[protocol]:
                    # print(results[protocol])
                    flat = [g for l in results[protocol][params]["gaps"] for g in l]
                    if flat:
                        # avg_island_len = sum(flat)/len(flat)
                        # print(protocol, params)
                        e_size = get_e_size(flat, args.L, args.nr_exp)
                    # else:
                        # avg_island_len = 0
                    res = [
                        round(100*results[protocol][params]["m"]/results[protocol][params]["mp"], 1),
                        100*results[protocol][params]["sc"]/(args.L*args.nr_exp),
                        100*results[protocol][params]["mc"]/(args.L*args.nr_exp),
                        e_size
                    ]
                    print(protocol, " & ", params, " & ", " & ".join([str(round(r, 1)) for r in res]), " & ", mut_freq)

    # print(results)

    # # random mutation positions
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
    parser.add_argument('--L', type=int, default=10000, help='Length of simulated sequences')
    parser.add_argument('--nr_exp', type=int, default=1000, help='Number of simulated experiments')
    parser.add_argument('--experiment_type', type=str, default="all", help='experiment type choose between "all", "controlled", "specified" or "only_subs"')
    parser.add_argument('--subs_freq', type=float, default=0.33, help='substitution frequency among all mutations for --experiment_type "specified"; rest split evenly in insertions and deletions')
    parser.add_argument('--mut_freqs', nargs='+', type=float, default=[0.01, 0.05, 0.10], help='mutation frequencies [0,1]')
    parser.add_argument('--k_size', type=int, default=30, help='k-mer/strobemer length')
    parser.add_argument('--w', type=int, default=1, help='number of hashes used in a sliding window for thinning (w=1 means no thinning)')
    #parser.add_argument('--orders', type=list, default=[4, ], help='List with orders of strobes to be analzyed')
    parser.add_argument('--orders', nargs='+', type=int, default=[2, ], help='List with orders of strobes to be analzyed')
    parser.add_argument('--w_low', type=int, default=25, help='minimum window offset to the previous window (wMin > 0)')
    parser.add_argument('--w_high', type=int, default=50, help='maximum window offset to the previous window (wMin <= wMax)')
    parser.add_argument('--strobe_fractions', nargs='+', type=float, default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], help='Fraction of sampled strobemers, rest kmers')
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

    methods = (
        "kmers", "spaced_kmers_dense", "spaced_kmers_sparse",
        "minstrobes", "randstrobes", "hybridstrobes", "altstrobes",
        "mixedminstrobes", "mixedrandstrobes", "mixedhybridstrobes", "mixedaltstrobes")

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
