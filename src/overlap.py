#!/usr/bin/env python3.9
# -*- coding: utf-8 -*
import os
import sys
import argparse
import random
import numpy as np
import collections
import seaborn as sns
import pandas as pd

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

from fractions import Fraction
from matplotlib import pyplot
from itertools import chain
from typing import Iterator
from typing import Union

from modules import indexing_Maier_altstrobes as indexing


def compute_overlap(seeding_technique: str, strobemers: list, nr_exp: int,
                    m: int, strobe_fraction: Union[str, int], n: int, outfile: str):
    """
    Compute weighted overlap between neighboring seeds

    :param seeding_technique: string with strobemer seeding technique
    :param strobemers: list with strobemer positions
    :param nr_exp: number of strobemers to be analyzed
    :param m: length of each individual strobe
    :param strobe_fraction: fraction of sampled strobemers, rest kmers
    :param n: number of neighboring seeds to be analyzed
    :param outfile: string with location and name of output file
    :returns: outputs overlap results in outfile
    """

    overlap = 0  # reset overlap
    for pos in range(nr_exp):
        strobe1 = set(range(strobemers[pos][0], strobemers[pos][0]+m)) | set(range(strobemers[pos][1], strobemers[pos][1]+m))
        for i in range(1, n+1):
            strobe2 = set(range(strobemers[pos+i][0], strobemers[pos+i][0]+m)) | set(range(strobemers[pos+i][1], strobemers[pos+i][1]+m))
            overlap += max(0, len(strobe1 & strobe2))/i
    overlap_per_exp = overlap/nr_exp
    with open(outfile, "a") as myfile:
        myfile.write(
            seeding_technique + " & " +
            str(strobe_fraction) + " & " +
            str(overlap_per_exp) + " & " +
            str(n) + "\n"
        )


def altstrobes(nr_exp: int, k: int, strobe_lengths: list, w_low: int,
               w_high: int, w: int, order: int, nr_neighboring: list, outfile: str):
    """
    Computes the weighted overlap of neighboring seeds from simulated sequences

    :param nr_exp: number of strobemers to be analyzed
    :param k: length of all strobes combined
    :param strobe_lengths: list with all altstrobe combinations (k1+k2=k) to be analyzed
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param nr_neighboring: number of neighboring seeds to be included in the analysis
    :param outfile: string with location and name of output file
    """

    for n in nr_neighboring:
        #assert n < k, "Number of neighboring seeds has to be smaller than k in this implementation"
        overlap_kmers = 0
        for i in range(1, n+1):
            overlap_kmers += max((k-i)/i, 0)
        with open("altstrobes_overlap.txt", "a") as myfile:
            myfile.write(
                "kmers & " +
                str(k) + " & " +
                str(overlap_kmers) + " & " +
                str(n) + "\n"
            )

        seq = "".join([random.choice("ACGT") for i in range(nr_exp + 2*order*w_high+k)])
        for k1, k2 in strobe_lengths:
            assert k1 + k2 == k, "[Warning] Lengths of strobes ({0}) is not adding up to k ({1})".format(k1+k2, k)
            overlap = 0
            strobemers = list(indexing.altstrobes_generalized(seq, k1, k2, w_low, w_high, w, order).keys())
            for pos in range(nr_exp):
                if strobemers[pos][0] != strobemers[pos][1]:  # (k1, k2)
                    strobe1 = set(range(strobemers[pos][0], strobemers[pos][0]+k1)) | set(range(strobemers[pos][2], strobemers[pos][2]+k2))
                else:  # (k2, k1)
                    strobe1 = set(range(strobemers[pos][0], strobemers[pos][0]+k2)) | set(range(strobemers[pos][2], strobemers[pos][2]+k1))

                for i in range(1, n+1):
                    if strobemers[pos+i][0] != strobemers[pos+i][1]:  # (k1, k2)
                        strobe2 = set(range(strobemers[pos+i][0], strobemers[pos+i][0]+k1)) | set(range(strobemers[pos+i][2], strobemers[pos+i][2]+k2))
                    else:  # (k2, k1)
                        strobe2 = set(range(strobemers[pos+i][0], strobemers[pos+i][0]+k2)) | set(range(strobemers[pos+i][2], strobemers[pos+i][2]+k1))

                    overlap += max(0, len(strobe1 & strobe2))/i
            overlap_per_exp = overlap/nr_exp
            if k1 == k2:
                with open("altstrobes_overlap.txt", "a") as myfile:
                    myfile.write(
                        "randstrobes & " +
                        str(k2) + " & " +
                        str(overlap_per_exp) + " & " +
                        str(n) + "\n"
                    )
            else:
                with open("altstrobes_overlap.txt", "a") as myfile:
                    myfile.write(
                        "altstrobes & " +
                        str(k2) + " & " +
                        str(overlap_per_exp) + " & " +
                        str(n) + "\n"
                    )


def mixedstrobes(nr_exp: int, k: int, w_low: int, w_high: int, w: int,
                 order: int, strobe_fractions: list, nr_neighboring: list,
                 outfile: str):
    """
    Computes the weighted overlap of neighboring seeds from simulated sequences

    :param k: length of all strobes combined
    :param w_low: minimum window offset to the previous window (wMin > 0)
    :param w_high: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fractions: list with all strobe_fractions (strobemer/k-mer) to be analyzed
    :param nr_neighboring: number of neighboring seeds to be included in the analysis
    :param outfile: string with location and name of output file
    """

    for n in nr_neighboring:
        #assert n < k, "Number of neighboring seeds has to be smaller than k in this implementation"
        overlap_kmers = 0
        for i in range(1, n+1):
            overlap_kmers += max((k-i)/i, 0)
        with open(outfile, "a") as myfile:
            myfile.write("kmers & " + str(0) + " & " + str(overlap_kmers) + " & " + str(n) + "\n")

        seq = "".join([random.choice("ACGT") for i in range(nr_exp + 2*order*w_high+k)])
        m = int(k//order)
        for s in strobe_fractions:

            # analyze mixedrandstrobes
            mixedstrobes = list(indexing.mixedrandstrobes(seq, k, w_low, w_high, w, order, s).keys())
            compute_overlap("randstrobes", mixedstrobes, nr_exp, m, s, n, outfile)

            # analyze mixedhybridstrobes
            mixedstrobes = list(indexing.mixedhybridstrobes(seq, k, w_low, w_high, w, order, s).keys())
            compute_overlap("hybridstrobes", mixedstrobes, nr_exp, m, s, n, outfile)

            # analyze mixedminstrobes
            mixedstrobes = list(indexing.mixedminstrobes(seq, k, w_low, w_high, w, order, s).keys())
            compute_overlap("minstrobes", mixedstrobes, nr_exp, m, s, n, outfile)


if __name__ == '__main__':
    k = 30
    strobe_lengths = [
        (15,15), (16,14), (17,13), (18,12), (19,11), (20,10), (21,9),
        (22,8), (23,7), (24,6), (25,5), (26,4), (27,3), (28,2), (29,1)
    ]
    strobe_fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    nr_exp = 1
    order = 2
    w_low = 25  # minimum strobe offset
    w_high = 50  # maximum strobe offset
    w = 1  # w = 1 means no thinning
    nr_neighboring = [1, 5, 25, 80]

    altstrobes(nr_exp, k, strobe_lengths, w_low, w_high, w, order, nr_neighboring, outfile="altstrobes_overlap.txt")
    mixedstrobes(nr_exp, k, w_low, w_high, w, order, strobe_fractions, nr_neighboring, outfile="mixedstrobes_overlap.txt")
