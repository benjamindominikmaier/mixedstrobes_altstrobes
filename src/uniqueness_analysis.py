#!/usr/bin/env python3.9
# -*- coding: utf-8 -*
import os
import sys
import argparse
import random

from collections import defaultdict
from modules import indexing_Maier_altstrobes as indexing, help_functions


def print_stats(acc: str, datastructure: str, all_mers: dict, k_size: int, total_mers: int) -> None:
    """
    Print uniqueness and ehits based on generated seeds for a given sequence

    :param acc: string containing the accession number
    :param datastructure: string with seeding technique
    :param all_mers: dictionary with all hash values as values
    :param k_size: length of the kmer/strobes
    :param total_mers: number of seeds for a given sequence and k-size
    """
    abundances = list(all_mers.values())
    del all_mers
    unique = abundances.count(1)
    percent_unique = round(100*unique/total_mers, 1)
    mean = sum(abundances)/len(abundances)
    ab_sorted = sorted(abundances)
    lower_75 = ab_sorted[1*len(abundances)//4]
    median = ab_sorted[len(abundances)//2]
    upper_75 = ab_sorted[3*len(abundances)//4]
    ehits = round(sum(i*i for i in abundances)/sum(abundances),2)

    data = ",".join([
        str(d)
        for d in [datastructure, k_size, acc, mean, median,
                  lower_75, upper_75, percent_unique, ehits]
    ])
    print(data)


def compute_uniqueness(args, acc: str, seq: str, k_size: int, w_low: int,
                       w_high: int, w: int, total_mers: int, order: int,
                       fraction: float=1, k_size2: int=0) -> None:
    """
    Compute uniqueness and e-hits of seeds for a given sequence and seeding technique

    :param args: user-specified arguments from argparser
    :param acc: string containing the accession number
    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmer/strobes
    :param total_mers: number of seeds for a given sequence and k-size
    :param order: number of substrings/strobes
    :param fraction: fraction of sampled strobemers, rest kmers
    """

    all_mers = defaultdict(int)
    if args.kmers:
        datastructure = "kmers"
        for i in range(len(seq) - k_size + 1):
            all_mers[hash(seq[i:i+k_size])] += 1

    if args.spaced_dense:
        datastructure = "spaced_dense"
        span_size = k_size+k_size//2
        positions = set(random.sample(range(1, span_size - 1), k_size-2))
        # asserts first and last position is sampled so that we have a spaced kmer of length span size
        positions.add(0)
        positions.add(span_size - 1)
        for s in indexing.spaced_kmers_iter(seq, k_size, span_size, positions):
            all_mers[s] += 1

    if args.spaced_sparse:
        datastructure = "spaced_sparse"
        span_size = 3*k_size
        positions = set(random.sample(range(1, span_size - 1), k_size-2))
        # asserts first and last position is sampled so that we have a spaced kmer of length span size
        positions.add(0)
        positions.add(span_size - 1)
        for s in indexing.spaced_kmers_iter(seq, k_size, span_size, positions):
            all_mers[s] += 1

    elif args.minstrobes:
        datastructure = "minstrobes" + str(order)
        for _, s in indexing.minstrobes_iter(seq, k_size, w_low, w_high, w, order, buffer_size=10000000):
            all_mers[s] += 1

    elif args.mixedminstrobes:
        datastructure = "mixedminstrobes" + str(order)
        for _, s in indexing.mixedminstrobes_iter(seq, k_size, w_low, w_high, w, order, fraction, buffer_size=10000000):
            all_mers[s] += 1

    elif args.randstrobes:
        datastructure = "randstrobes" + str(order)
        for _, s in indexing.randstrobes_iter(seq, k_size, w_low, w_high, w, order, buffer_size=10000000):
            all_mers[s] += 1

    elif args.mixedrandstrobes:
        datastructure = "mixedrandstrobes" + str(order)
        for _, s in indexing.mixedrandstrobes_iter(seq, k_size, w_low, w_high, w, order, fraction, buffer_size=10000000):
            all_mers[s] += 1

    elif args.hybridstrobes:
        datastructure = "hybridstrobes" + str(order)
        for _, s in indexing.hybridstrobes_iter(seq, k_size, w_low, w_high, w, order, buffer_size=10000000):
            all_mers[s] += 1

    elif args.mixedhybridstrobes:
        datastructure = "mixedhybridstrobes" + str(order)
        for _, s in indexing.mixedhybridstrobes_iter(seq, k_size, w_low, w_high, w, order, fraction, buffer_size=10000000):
            all_mers[s] += 1

    elif args.altstrobes:
        datastructure = "altstrobes" + str(order)
        m_size1 = int(k_size/3)
        m_size2 = int(2*k_size/3)

        k_size = (m_size1, m_size2)
        for _, s in indexing.altstrobes_iter(seq, m_size1, m_size2, w_low, w_high, w, order, buffer_size=10000000):
            all_mers[s] += 1

    elif args.mixedaltstrobes:
        datastructure = "mixedaltstrobes" + str(order)
        m_size1 = int(k_size/3)
        m_size2 = int(2*k_size/3)
        k_size = (m_size1, m_size2)
        for _, s in indexing.mixedaltstrobes_iter(seq, m_size1, m_size2, w_low, w_high, w, order, fraction, buffer_size=10000000):
            all_mers[s] += 1

    print_stats(acc, datastructure, all_mers, k_size, total_mers)


def altstrobes_generalized_uniqueness(args):
    """
    Helper to compute the uniqueness of seeds yielded from given or simulated sequences
    using generalized altstrobes where k2 != 2*k1
    """
    if args.fasta:
        genome = {acc: seq for (acc, (seq, _)) in help_functions.readfq(open(args.fasta, 'r'))}
        for acc, seq in genome.items():
            acc = acc.split()[0]
            # print(acc)
            genome[acc] = seq.replace("N", "")  # remove Ns

        # print(len(genome), sum([len(v) for k, v in genome.items()]))
        print("datastructure,k,acc,mean,median,lower_75,upper_75,\%-unique,ehits")
        total_mers = {}
        for acc, seq in genome.items():
            total_mers[acc] = len(seq) - args.altstrobes_generalized + 1
            # if acc == "chr1" or acc == "chr2" or acc == "chr3":  # or acc == "chr4" or acc == "chr5":
            #if acc == "NC_000001.11" or acc == "NC_000002.12" or acc == "NC_000003.12":  # or acc == "chr4" or acc == "chr5":
            #    chr_dict = {"NC_000001.11": "chr1", "NC_000002.12": "chr2", "NC_000003.12": "chr3"}
            if acc == "NC_000021.9":  # chr21
            # if acc == "NC_000001.11":  # chr1
                for k1 in range(1, args.altstrobes_generalized+1):
                    k_size = (k1, args.altstrobes_generalized-k1)
                    datastructure = "altstrobes" + str(order) + str(k_size)
                    for _, s in indexing.altstrobes_iter(seq, k_size[0], k_size[1], args.w_low, args.w_high, args.w, args.order, buffer_size=10000000):
                        all_mers[s] += 1
                    print_stats(acc, datastructure, all_mers, k_size, total_mers[acc])
    else:
        print("Simulating sequence of length: ", args.L)
        for exp_id in range(args.nr_exp):
            seq = "".join([random.choice("ACGT") for i in range(args.L)])
            for k1 in range(1, args.altstrobes_generalized+1):
                k_size = (k1, args.altstrobes_generalized-k1)
                datastructure = "altstrobes" + str(order) + str(k_size)
                for _, s in indexing.altstrobes_iter(seq, k_size[0], k_size[1], args.w_low, args.w_high, args.w, args.order, buffer_size=10000000):
                    all_mers[s] += 1
                print_stats(acc, datastructure, all_mers, k_size, total_mers)


def uniqueness(args):
    """
    Helper to compute the uniqueness of seeds yielded from given or simulated sequences
    using strobemers or k-mers
    """
    if args.fasta:
        genome = {acc: seq for (acc, (seq, _)) in help_functions.readfq(open(args.fasta, 'r'))}
        for acc, seq in genome.items():
            acc = acc.split()[0]
            # print(acc)
            genome[acc] = seq.replace("N", "")  # remove Ns

        # print(len(genome), sum([len(v) for k, v in genome.items()]))
        print("datastructure,k,acc,mean,median,lower_75,upper_75,\%-unique,ehits")

        total_mers = {}
        for acc, seq in genome.items():
            # print(acc)
            for k_size in args.k_sizes:  # [18,24,30,36]:
                total_mers[acc] = len(seq) - k_size + 1
                # if acc == "chr1" or acc == "chr2" or acc == "chr3":  # or acc == "chr4" or acc == "chr5":
                #if acc == "NC_000001.11" or acc == "NC_000002.12" or acc == "NC_000003.12":  # or acc == "chr4" or acc == "chr5":
                #    chr_dict = {"NC_000001.11": "chr1", "NC_000002.12": "chr2", "NC_000003.12": "chr3"}
                if acc == "NC_000021.9":  # chr21
                # if acc == "NC_000001.11":  # chr1
                    compute_uniqueness(args, acc, seq, k_size, args.w_low, args.w_high, args.w, total_mers[acc], args.order, args.strobe_fraction)
    else:
        print("Simulating sequence of length: ", args.L)
        for exp_id in range(args.nr_exp):
            seq = "".join([random.choice("ACGT") for i in range(args.L)])
            for k_size in args.k_sizes:  # [18,24,30,36]:
                total_mers = len(seq) - k_size + 1
                compute_uniqueness(args, "", seq, k_size, args.w_low, args.w_high, args.w, total_mers, args.order, args.strobe_fraction)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    parser.add_argument('--kmers',  action="store_true", help='Seeding k-mers')
    parser.add_argument('--minstrobes', action="store_true", help='Seeding minstrobes')
    parser.add_argument('--mixedminstrobes',  action="store_true", help='Seeding mixedminstrobes (minstrobes/k-mers) with a user-defined --strobe_fraction')
    parser.add_argument('--randstrobes',  action="store_true", help='Seeding randstrobes')
    parser.add_argument('--mixedrandstrobes',  action="store_true", help='Seeding mixedrandstrobes (randstrobes/k-mers) with a user-defined --strobe_fraction')
    parser.add_argument('--hybridstrobes',  action="store_true", help='Seeding hybridstrobes')
    parser.add_argument('--mixedhybridstrobes',  action="store_true", help='Seeding mixedhybridstrobes (hybridstrobes/k-mers) with a user-defined --strobe_fraction')
    parser.add_argument('--spaced_dense',  action="store_true", help='Seeding spaced k-mers (dense)')
    parser.add_argument('--spaced_sparse',  action="store_true", help='Seeding spaced k-mers (sparse)')
    parser.add_argument('--altstrobes',  action="store_true", help='Kmer size')
    parser.add_argument('--mixedaltstrobes', action="store_true", help='Kmer size')
    parser.add_argument('--order', type=int, default=2, help='Order on strobes')
    parser.add_argument('--strobe_fraction', type=float, default=1, help='Fraction of sampled strobemers, rest kmers')
    parser.add_argument('--k_sizes', '--nargs-int-type', nargs='+', type=int, default=[18, 24, 30, 36], help="List with strobe lengths to be analyzed")
    parser.add_argument('--w', type=int, default=1, help="number of hashes used in a sliding window for thinning (w=1 means no thinning)")
    parser.add_argument('--w_low', type=int, default=25, help="minimum window offset to the previous window (wMin > 0)")
    parser.add_argument('--w_high', type=int, default=50, help="maximum window offset to the previous window (wMin <= wMax)")
    parser.add_argument('--altstrobes_generalized', type=int, default=0, help='Choose k-size to seed altstrobes with strobe combinations from (1,k-1) to (k-1,1)')
    parser.add_argument('--L', type=int, default=1000, help='Length of simulated sequences (only if no fasta file is provided)')
    parser.add_argument('--nr_exp', type=int, default=1, help='Number of simulated experiments (only if no fasta file is provided)')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    if args.altstrobes_generalized > 0:
        altstrobes_generalized_uniqueness(args)
    else:
        uniqueness(args)
