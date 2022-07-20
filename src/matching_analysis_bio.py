#! /usr/bin/env python

from __future__ import print_function
import os
import sys
import argparse
import copy
import operator
import math

# import errno
# from time import time
# import re

# import random
# import parasail
# import pysam

from collections import defaultdict, deque
from sys import stdout
from array import array
from itertools import zip_longest
from typing import Iterator
from fractions import Fraction

from modules import help_functions
from genome_mapping_metrics import *

BITS = sys.hash_info.width
MAX = sys.maxsize
MAX_HASH_VALUE = int((2**BITS)/2) - 1


def argmin(array: list) -> tuple:
    """
    Find the value of x which minimizes f(x) over the set of candidates for x

    :param array: a list to minimize
    :returns: a tuple with the index position and the value of the lowest element
    """
    min_index = array.index(min(array))
    min_val = array[min_index]
    return min_index, min_val


def reverse_complement(seq: str) -> str:
    """
    Compute the reverse complement of a given nucleotide sequence

    :param seq: a string with a nucleotide sequence
    :returns: a string with the reverse complement of the given sequence
    """
    rev_nuc = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
        'N': 'N', 'X': 'X', 'n': 'n',
        'Y': 'R', 'R': 'Y', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V',
        'V': 'B', 'H': 'D', 'D': 'H',
        'y': 'r', 'r': 'y', 'k': 'm', 'm': 'k', 's': 's', 'w': 'w', 'b': 'v',
        'v': 'b', 'h': 'd', 'd': 'h'
    }
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(seq)])
    return(rev_comp)


def thinner(hash_list: list, w: int) -> list:
    """
    Thins out kmers/strobemers using a sliding window approach

    :param hash_list: a list with hash values
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: a list with tuples (pos in original list, minimim hash value) for each window of w hashes
    """
    window_hashes = deque(hash_list[:w])
    min_index, curr_min_hash = argmin(window_hashes)
    thinned_hash_list = [(min_index, curr_min_hash)]

    for i in range(w, len(hash_list) + w-1):
        if i >= len(hash_list):
            new_hash = MAX
        else:
            new_hash = hash_list[i]

        # updating window
        discarded_hash = window_hashes.popleft()
        window_hashes.append(new_hash)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min_hash == discarded_hash:
            min_index, curr_min_hash = argmin(window_hashes)
            thinned_hash_list.append((min_index + i + 1 - w, curr_min_hash))

        # previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_hash < curr_min_hash:
            curr_min_hash = new_hash
            thinned_hash_list.append((i, curr_min_hash))

    return thinned_hash_list


def update_queue(q: list, curr_min: int, min_index: int, new_hash: int, i: int,
                 start_offset: int, end_offset: int) -> tuple:
    """
    Updates windows

    :param q: a list with strobe_windows
    :param curr_min: an integer with the current minimum value
    :param min_index: an integer with the index position of the minimum value
    :param new_hash: an integer with the new hash value
    :param i: an integer with the position of the first strobe
    :param start_offset: minimum window offset
    :param end_offset: maximum window offset
    :returns: a tuple with the index position of the minimum value and the value
    """
    old_h = q.popleft()
    q.append(new_hash)

    # we have discarded previous windows minimizer, look for new minimizer brute force
    if curr_min == old_h:
        min_index, curr_min = argmin(q)
        min_index = i + start_offset + min_index

    # Previous minimizer still in window, we only need to compare with the recently added kmer
    elif new_hash < curr_min:
        curr_min = new_hash
        min_index = i + end_offset

    return min_index, curr_min


def seq_to_hybridstrobes_iter(seq: str, k_size: int, w_min, w_max, prime: int,
                              w: int, order: int) -> Iterator[tuple]:
    """
    Generator for creating hybridstrobes of any orders

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param w_min: minimum window offset to the previous window (wMin > 0)
    :param w_max: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: an iterator for creating hybridstrobes
    """
    if len(seq) < 2*w_max:
        return (*(-1 for i in range(order)), None)

    hash_list = [
        hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)
        if "N" not in seq[i:i+k_size]]

    n_partition = 3
    w_p = (w_max - w_min) // n_partition

    tmp_index_dict = dict()
    for strobe_num in range(0, order-1):
        tmp_index = []
        for partition in range(0, n_partition):
            start = w_min + w_max*strobe_num + w_p*partition
            end = (
                w_max + w_max*strobe_num if partition + 1 == n_partition
                else w_min + w_max*strobe_num + w_p + w_p*partition
            )

            strobe_window = deque(hash_list[start: end])
            min_index, min_w = argmin(strobe_window)
            min_index = min_index + w_min + w_max*strobe_num + w_p*partition
            tmp_index.append(
                (
                    strobe_window,
                    min_w,
                    min_index,
                    start, end
                )
            )
        tmp_index_dict[strobe_num] = tmp_index

    for i in range(len(hash_list) - w_max*order):  # temporary iteration
        index_hash = hash_list[i]
        positions = [i, ]

        for strobe_num in range(0, order-1):
            tmp_index = []
            for window_numer, window in enumerate(tmp_index_dict[strobe_num]):
                # updating windows
                strobe_window, min_w, min_index, start, end = window
                new_w = hash_list[i + end]
                min_index, min_w = update_queue(
                    strobe_window, min_w, min_index, new_w, i, start, end
                )

                # update tmp_index_dict
                tmp_index_dict[strobe_num][window_numer] = (
                    strobe_window, min_w, min_index, start, end
                )
                tmp_index.append((min_index, min_w))

            next_i, next_index_hash = tmp_index[index_hash % n_partition]
            positions.append(next_i)
            index_hash = index_hash + (strobe_num+1) * (-1)**(strobe_num+1) * next_index_hash

        yield positions, index_hash


def altstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
               strobe_w_max_offset: int, w: int, order: int = 2,
               prime: int = 997) -> dict:
    """
    Strobemer seeding protocol to sample altstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :returns: a dictionary with positions along the string as keys and the altstrobes as value
    """

    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % (order+1) != 0:
        print("WARNING: kmer size {0} is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size1 = int(k_size/3)
    m_size2 = int(2*k_size/3)

    altstrobes = {tuple(index): h for index, h, min_values in seq_to_altstrobes_iter(
        seq, m_size1, m_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order
    )}
    return altstrobes

def altstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                     strobe_w_max_offset: int, w: int, order: int = 2,
                     buffer_size: int = 10000000) -> Iterator[tuple]:
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in altstrobes(
                substring, k_size, strobe_w_min_offset, strobe_w_max_offset,
                w, order=order).items():

            yield p, m

def seq_to_mixedaltstrobes_iter(seq: str, k_size1: int, k_size2: int, strobe_w_min_offset: int,
                                strobe_w_max_offset: int, prime: int, w: int,
                                order: int, denominator: int, numerator: int) -> Iterator[tuple]:
    """
    Iterator for creating of mixedaltstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: an iterator for creating mixedaltstrobes
    """
    hash_seq_list1 = [(i, hash(seq[i:i+k_size1])) for i in range(len(seq) - k_size1 + 1)]
    hash_seq_list2 = [(i, hash(seq[i:i+k_size2])) for i in range(len(seq) - k_size2 + 1)]
    strobe_w_max_offset2 = strobe_w_max_offset - k_size1 # 40 vs 30
    strobe_w_min_offset2 = strobe_w_min_offset - k_size1 # 21 vs 11

    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned1 = thinner([h for i, h in hash_seq_list1], w)
        hash_seq_list_thinned2 = thinner([h for i, h in hash_seq_list2], w)
    else:
        hash_seq_list_thinned1 = hash_seq_list1
        hash_seq_list_thinned2 = hash_seq_list2

    for index_position, (p1, hash_m1) in enumerate(hash_seq_list_thinned1):  # [:-k_size]:
        if p1 >= len(hash_seq_list2) - (order-1)*k_size2:
            break
        # hash_m1 = hash_seq_list[p]

        if hash_m1 % denominator < numerator:  # pick altstrobe
            if hash_m1 % 2 == 0: # first x, then 2x (e.g. 10-20)
                # print("10-20")
                if p1 + (order-1) * strobe_w_max_offset2 <= len(hash_seq_list2):
                    windows = list()
                    for window_order in range(1, order):
                        start = p1 + strobe_w_min_offset2 + (window_order-1) * strobe_w_max_offset2
                        end = min(p1 + window_order * strobe_w_max_offset2, len(hash_seq_list2))
                        windows.append((start, end))

                else:
                    windows = list()
                    for window_order in range(1, order):
                        start = (max(
                            p1+window_order*k_size2,
                            len(hash_seq_list2) + strobe_w_min_offset2 - (order - window_order) * strobe_w_max_offset2
                            )
                        )

                        end = min(p1 + window_order * strobe_w_max_offset2, len(hash_seq_list2))
                        windows.append((start, end))

                index = [p1, ] # -p1
                min_values = []
                min_hash_val = hash_m1
                for index_order in range(1, order):
                    min_index, min_value = argmin([
                        (min_hash_val + hash_seq_list2[i][1]) % prime
                        for i in range(*windows[index_order-1])
                    ])

                    min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list2[windows[index_order-1][0] + min_index][1]
                    index.append(min_index+windows[index_order-1][0])
                    index.append(min_index+windows[index_order-1][0]+k_size1)
                    min_values.append(min_value)

            else: # first 2x, then x (e.g. 20-10)
                # print("20-10")
                if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list1):
                    windows = list()
                    for window_order in range(1, order):
                        start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                        end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list1))
                        windows.append((start, end))

                else:
                    windows = list()
                    for window_order in range(1, order):
                        start = (max(
                            p1+window_order*k_size1,
                            len(hash_seq_list1) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                            )
                        )

                        end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list1))
                        windows.append((start, end))

                index = [p1, p1+k_size1]
                min_values = []
                min_hash_val = hash_seq_list2[index_position][1]
                for index_order in range(1, order):
                    min_index, min_value = argmin([
                        (min_hash_val + hash_seq_list1[i][1]) % prime
                        for i in range(*windows[index_order-1])
                    ])

                    min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list1[windows[index_order-1][0] + min_index][1]
                    index.append(min_index+windows[index_order-1][0]) # -

            yield index, min_hash_val

        else:  # pick k-mer
            index = tuple(p1 + (strobe_num) * k_size1 for strobe_num in range(order+1))
            l = int(k_size1+k_size2/2*order)
            yield index, hash(seq[p1: p1+l])


def mixedaltstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int, order: int = 2,
                    strobe_fraction: float = 0.5) -> dict:
    """
    Strobemer seeding protocol to sample mixedaltstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :returns: a dictionary with positions along the string as keys and the altstrobes as value
    """
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % (order+1) != 0:
        print("WARNING: kmer size is not evenly divisible with {0}, will use {1} as kmer size: ".format(order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size1 = int(k_size/3)
    m_size2 = int(2*k_size/3)

    mixedaltstrobes = {tuple(index): h for index, h in seq_to_mixedaltstrobes_iter(
        seq, m_size1, m_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator
    )}

    return mixedaltstrobes


def seq_to_randstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                            strobe_w_max_offset: int, prime: int, w: int,
                            order: int) -> Iterator[tuple]:
    """
    Iterator for creation of randstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: an iterator for creating randstrobes
    """
    hash_seq_list = [
        (i, hash(seq[i: i+k_size])) for i in range(len(seq) - k_size + 1)
        if "N" not in seq[i: i+k_size]
    ]
    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 + (order-1) * strobe_w_max_offset > len(hash_seq_list):
            break
        # hash_m1 = hash_seq_list[p]

        else:
            windows = list()
            for window_order in range(1, order):
                start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                windows.append((start, end))

        positions = [p1, ]
        min_values = []
        min_hash_val = hash_m1
        for index_order in range(1, order):
            min_index, min_value = argmin([
                (min_hash_val + hash_seq_list[i][1]) % prime
                for i in range(*windows[index_order-1])
            ])

            # min_values.append(min_value)
            min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list[windows[index_order-1][0] + min_index][1]
            positions.append(min_index+windows[index_order-1][0])

        yield positions, min_hash_val


def seq_to_minstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                           strobe_w_max_offset: int, prime: int, w: int,
                           order: int) -> Iterator[tuple]:
    """
    Generator for creating minstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: an iterator for creating minstrobes
    """
    hash_seq_list = [
        (i, hash(seq[i: i+k_size])) for i in range(len(seq) - k_size + 1)
        if "N" not in seq[i:i+k_size]
    ]

    # produce a subset of positions, still with samme index as in full sequence
    strobes = deque(thinner([h for i, h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset))
    strobes_dict = {strobe_num: copy.deepcopy(strobes) for strobe_num in range(1, order)}

    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-m_size]:
        if p1 >= len(hash_seq_list) + k_size - k_size*order:
            break

        positions = [p1, ]
        hash_value = hash_m1

        for strobe_num in range(1, order):
            if p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset < len(seq):
                while strobes_dict[strobe_num][0][0] < min(p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset, len(hash_seq_list)-1):
                    l = strobes_dict[strobe_num].popleft()
            p, hash_val = strobes_dict[strobe_num][0]
            positions.append(p)
            hash_value += strobe_num * (-1)**strobe_num * hash_val

        yield positions, hash_value


def seq_to_mixedhybridstrobes_iter(seq: str, k_size: int, w_min: int,
                                   w_max: int, prime: int, w: int, order: int,
                                   denominator: int,
                                   numerator: int) -> Iterator[tuple]:
    """
    Generator for creating mixed hybridstrobes of any orders

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param w_min: minimum window offset to the previous window (wMin > 0)
    :param w_max: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: an iterator for creating mixedhybridstrobes
    """
    hash_list = [
        hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)
        if "N" not in seq[i:i+k_size]
    ]

    n_partition = 3
    w_p = (w_max - w_min) // n_partition

    tmp_index_dict = dict()
    for strobe_num in range(0, order-1):
        tmp_index = []
        for partition in range(0, n_partition):
            start = w_min + w_max*strobe_num + w_p*partition
            end = (
                w_max + w_max*strobe_num if partition + 1 == n_partition
                else w_min + w_max*strobe_num + w_p + w_p*partition
            )

            strobe_window = deque(hash_list[start: end])
            min_index, min_w = argmin(strobe_window)
            min_index = min_index + w_min + w_max*strobe_num + w_p*partition
            tmp_index.append(
                (
                    strobe_window,
                    min_w,
                    min_index,
                    start, end
                )
            )
        tmp_index_dict[strobe_num] = tmp_index

    for i in range(len(hash_list) - w_max*order):  # temporary iteration
        index_hash = hash_list[i]
        positions = [i, ]

        for strobe_num in range(0, order-1):
            tmp_index = []
            for window_numer, window in enumerate(tmp_index_dict[strobe_num]):
                # updating windows
                strobe_window, min_w, min_index, start, end = window
                new_w = hash_list[i + end]
                min_index, min_w = update_queue(
                    strobe_window, min_w, min_index, new_w, i, start, end
                )

                # update tmp_index_dict
                tmp_index_dict[strobe_num][window_numer] = (
                    strobe_window, min_w, min_index, start, end
                )
                tmp_index.append((min_index, min_w))

            next_i, next_index_hash = tmp_index[index_hash % n_partition]
            positions.append(next_i)
            index_hash = index_hash + (strobe_num+1) * (-1)**(strobe_num+1) * next_index_hash

        # decide whether kmers should be sampled instead of mixedstrobes
        if hash_list[i] % denominator >= numerator:
            positions = [i + strobe * k_size for strobe in range(order)]
            index_hash = hash(seq[i:i+k_size * order])

        yield positions, index_hash


def seq_to_mixedminstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                                strobe_w_max_offset: int, prime: int, w: int,
                                order: int, denominator: int,
                                numerator: int) -> Iterator[tuple]:
    """
    Generator for creating mixedminstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: a tuple with positions as first element and hash_value as second element.
    """
    hash_seq_list = [
        (i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)
        if "N" not in seq[i:i+k_size]
    ]

    # produce a subset of positions, still with samme index as in full sequence
    strobes = deque(thinner([h for i, h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset))
    strobes_dict = {strobe_num: copy.deepcopy(strobes) for strobe_num in range(1, order)}

    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    # Decision based on hash values
    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 >= len(hash_seq_list) + k_size - order*k_size:
            break
        if hash_m1 % denominator < numerator:  # pick minstrobe
            positions = [p1, ]
            hash_value = hash_m1

            for strobe_num in range(1, order):
                if p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset < len(seq):
                    while strobes_dict[strobe_num][0][0] < min(p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset, len(hash_seq_list)-1):
                        l = strobes_dict[strobe_num].popleft()
                p, hash_val = strobes_dict[strobe_num][0]
                positions.append(p)
                hash_value += strobe_num * (-1)**strobe_num * hash_val
            yield positions, hash_value

        else:  # pick k-mer
            positions = tuple(p1 + (strobe_num) * k_size for strobe_num in range(order))
            yield positions, hash(seq[p1:p1+k_size*order])


def seq_to_mixedrandstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                                 strobe_w_max_offset: int, prime: int, w: int,
                                 order: int, denominator: int,
                                 numerator: int) -> Iterator[tuple]:
    """
    Iterator for creating of mixedrandstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: an iterator for creating mixedrandstrobes
    """
    hash_seq_list = [
        (i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)
        if "N" not in seq[i:i+k_size]
    ]

    # thinning
    if w > 1:
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 + (order-1) * strobe_w_max_offset > len(hash_seq_list):
            break
        # hash_m1 = hash_seq_list[p]

        if hash_m1 % denominator < numerator:  # pick randstrobe
            windows = list()
            for window_order in range(1, order):
                start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                windows.append((start, end))


            positions = [p1, ]
            min_hash_val = hash_m1
            for index_order in range(1, order):
                min_index, min_value = argmin([
                    (min_hash_val + hash_seq_list[i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])

                min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list[windows[index_order-1][0] + min_index][1]
                positions.append(min_index+windows[index_order-1][0])
            yield positions, min_hash_val

        else:  # pick k-mer
            positions = tuple(p1 + (strobe_num) * k_size for strobe_num in range(order))
            yield positions, hash(seq[p1: p1+k_size*order])


def seq_to_kmer_iter(seq: str, k_size: int, w: int) -> Iterator[tuple]:
    """
    Iterator for creating kmers

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: an iterator for creating kmers
    """
    hash_seq_list = [
        (i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)
        if "N" not in seq[i:i+k_size]
    ]

    if w > 1:
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
        for p, h in hash_seq_list_thinned:
            yield ([p,], h)
    else:
        for p, h in hash_seq_list:
            yield ([p,], h)


def grouper(iterable, n: int, fillvalue=None) -> Iterator:
    """
    Collect data into fixed-length chunks or blocks without throwing excess items away

    :param iterable: sequence data type with values to combine
    :param n: length of chunks or blocks
    :param fillvalue: plugs the gap if there is no value to match to one of our iterable's elements
    :returns: a terminating iterator (zip_longest) with fixed-length chunks or blocks
    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def build_kmer_index(refs, k_size: int, w: int) -> tuple:
    """
    Build kmers from references

    :param refs: reference fasta file
    :param k_size: length of the kmer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: a tuple with an index-dictionary, a reference-to-accession-number dictionary and a counter of kmers created from references
    """
    idx = defaultdict(lambda: array("L"))
    ref_id_to_accession = {}
    cntr = 0
    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for positions, hash_val in seq_to_kmer_iter(seq, k_size, w):
            idx[hash_val].append(r_id)
            idx[hash_val].append(positions[0])
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} kmers created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def build_strobemer_index(method: str, refs, k_size: int,
                          strobe_w_min_offset: int, strobe_w_max_offset: int,
                          prime: int, w: int, order: int) -> tuple:
    """
    Build strobemers from references

    :param method: string which specifies whether minstrobes or randstrobes should be sampled
    :param refs: reference fasta file
    :param k_size: length of all strobes combined
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: rime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: a tuple with an index-dictionary, a reference-to-accession-number dictionary and a counter of strobemers created from references
    """
    idx = defaultdict(lambda: array("L"))
    ref_id_to_accession = {}
    cntr = 0
    method_function = "seq_to_" + method + "_iter"

    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for positions, hash_val in globals()[method_function](seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order):
            idx[hash_val].append(r_id)
            for pos in positions:
                idx[hash_val].append(pos)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} {1} created from references".format(cntr, method))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def build_mixedstrobemer_index(method: str, refs, k_size: int,
                               strobe_w_min_offset: int, strobe_w_max_offset: int,
                               prime: int, w: int, order: int,
                               denominator: int, numerator: int) -> tuple:
    """
    Build mixedstrobes (strobemer/kmer) from references

    :param method: string which specifies whether mixedminstrobes or mixedrandstrobes should be sampled
    :param refs: reference fasta file
    :param k_size: length of all strobes combined
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: a tuple with an index-dictionary, a reference-to-accession-number dictionary and a counter of mixedstrobemers created from references
    """
    idx = defaultdict(lambda: array("L"))
    ref_id_to_accession = {}
    cntr = 0
    method_function = "seq_to_" + method + "_iter"

    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for positions, hash_val in globals()[method_function](seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator):
            idx[hash_val].append(r_id)
            for pos in positions:
                idx[hash_val].append(pos)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} {1} created from references".format(cntr, method))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr


def build_hybridstrobe_index(refs, k_size: int, strobe_w_min_offset: int,
                             strobe_w_max_offset: int, prime: int, w: int,
                             order) -> tuple:
    """
    Build hybridstrobes from references

    :param refs: reference fasta file
    :param k_size: length of all strobes combined
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: a tuple with an index-dictionary, a reference-to-accession-number dictionary and a counter of hybridstrobes created from references
    """
    idx = defaultdict(lambda: array("L"))
    ref_id_to_accession = {}
    cntr = 0

    if w == 1:
        for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
            ref_id_to_accession[r_id] = ref_acc
            for positions, hash_val in seq_to_hybridstrobes_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order):
                idx[hash_val].append(r_id)
                for pos in positions:
                    idx[hash_val].append(pos)
                cntr += 1
                if cntr % 1000000 == 0:
                    print("{0} hybridstrobes created from references".format(cntr))
            # print(hash_val, r_id, pos)
    else:
        for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
            thinner_window = deque([hash(seq[i:i+k_size]) for i in range(w)])
            min_index, curr_min_hash = argmin(thinner_window)
            sampled_positions = set([min_index])
            ref_id_to_accession[r_id] = ref_acc
            info_buffer = deque([])
            for positions, hash_val in seq_to_hybridstrobes_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order):
                if positions[0] in sampled_positions:
                    idx[hash_val].append(r_id)
                    for pos in positions:
                        idx[hash_val].append(pos)
                    sampled_positions.remove(positions[0])
                    assert len(sampled_positions) == 0
                    cntr += 1

                if positions[0] < w:
                    info_buffer.append((positions, hash_val))
                    continue  # already in queue
                else:
                    # updating window
                    discarded_hash = thinner_window.popleft()
                    thinner_window.append(hash_val)
                    info_buffer.popleft()
                    info_buffer.append((positions, hash_val))

                    # we have discarded previous windows minimizer, look for new minimizer brute force
                    if curr_min_hash == discarded_hash:
                        min_index, curr_min_hash = argmin(thinner_window)
                        (positions_, hash_val_) = info_buffer[min_index]
                        idx[hash_val].append(r_id)
                        for pos in positions_:
                            idx[hash_val].append(pos)
                        cntr += 1

                    # Previous minimizer still in window, we only need to compare with the recently added kmer
                    elif hash_val < curr_min_hash:
                        curr_min_hash = hash_val
                        idx[hash_val].append(r_id)
                        for pos in positions:
                            idx[hash_val].append(pos)
                        cntr += 1

                if cntr % 1000000 == 0:
                    print("{0} hybridstrobes created from references, currently at position: {1}".format(cntr, positions[0]))

    return idx, ref_id_to_accession, cntr


def build_mixedhybridstrobe_index(refs, k_size: int,
                                  strobe_w_min_offset: int, strobe_w_max_offset: int,
                                  prime: int, w: int, order: int,
                                  denominator: int, numerator: int) -> tuple:
    """
    Build mixedhybridstrobes from references

    :param refs: reference fasta file
    :param k_size: length of all strobes combined
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns a tuple with an index-dictionary, a reference-to-accession-number dictionary and a counter of mixedhybridstrobes created from references
    """
    idx = defaultdict(lambda: array("L"))
    ref_id_to_accession = {}
    cntr = 0

    if w == 1:
        for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
            ref_id_to_accession[r_id] = ref_acc
            for positions, hash_val in seq_to_mixedhybridstrobes_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator):
                idx[hash_val].append(r_id)
                for pos in positions:
                    idx[hash_val].append(pos)
                cntr += 1
                if cntr % 1000000 == 0:
                    print("{0} hybridstrobes created from references".format(cntr))
            # print(hash_val, r_id, pos)
    else:
        for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
            thinner_window = deque([hash(seq[i:i+k_size]) for i in range(w)])
            min_index, curr_min_hash = argmin(thinner_window)
            sampled_positions = set([min_index])
            ref_id_to_accession[r_id] = ref_acc
            info_buffer = deque([])
            for positions, hash_val in seq_to_mixedhybridstrobes_iter(seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator):
                if positions[0] in sampled_positions:
                    idx[hash_val].append(r_id)
                    for pos in positions:
                        idx[hash_val].append(pos)
                    sampled_positions.remove(positions[0])
                    assert len(sampled_positions) == 0
                    cntr += 1

                if positions[0] < w:
                    info_buffer.append((positions, hash_val))
                    continue  # already in queue
                else:
                    # updating window
                    discarded_hash = thinner_window.popleft()
                    thinner_window.append(hash_val)
                    info_buffer.popleft()
                    info_buffer.append((positions, hash_val))

                    # we have discarded previous windows minimizer, look for new minimizer brute force
                    if curr_min_hash == discarded_hash:
                        min_index, curr_min_hash = argmin(thinner_window)
                        (positions_, hash_val_) = info_buffer[min_index]
                        idx[hash_val].append(r_id)
                        for pos in positions_:
                            idx[hash_val].append(pos)
                        cntr += 1

                    # Previous minimizer still in window, we only need to compare with the recently added kmer
                    elif hash_val < curr_min_hash:
                        curr_min_hash = hash_val
                        idx[hash_val].append(r_id)
                        for pos in positions:
                            idx[hash_val].append(pos)
                        cntr += 1

                if cntr % 1000000 == 0:
                    print("{0} hybridstrobes created from references, currently at position: {1}".format(cntr, positions[0]))

    return idx, ref_id_to_accession, cntr





def build_altstrobe_index(refs, k_size1: int, k_size2: int,
                          strobe_w_min_offset: int, strobe_w_max_offset: int,
                          prime: int, w: int, order: int) -> tuple:
    """
    Build altstrobes from references

    :param method: string which specifies whether minstrobes or randstrobes should be sampled
    :param refs: reference fasta file
    :param k_size: length of all strobes combined
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: rime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: a tuple with an index-dictionary, a reference-to-accession-number dictionary and a counter of strobemers created from references
    """
    idx = defaultdict(lambda: array("L"))
    ref_id_to_accession = {}
    cntr = 0

    k_size1 = 10
    k_size2 = 20

    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for positions, hash_val in seq_to_altstrobes_iter(seq, k_size1, k_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order):
            idx[hash_val].append(r_id)
            for pos in positions:
                idx[hash_val].append(pos)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} altstrobes created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr





def build_mixedaltstrobe_index(refs, k_size1: int, k_size2: int, strobe_w_min_offset: int, strobe_w_max_offset: int,
                          prime: int, w: int, order: int,
                          denominator: int, numerator: int) -> tuple:
    """
    Build mixedstrobes (strobemer/kmer) from references

    :param method: string which specifies whether mixedminstrobes or mixedrandstrobes should be sampled
    :param refs: reference fasta file
    :param k_size: length of all strobes combined
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: a tuple with an index-dictionary, a reference-to-accession-number dictionary and a counter of mixedstrobemers created from references
    """
    idx = defaultdict(lambda: array("L"))
    ref_id_to_accession = {}
    cntr = 0

    k_size1 = 10
    k_size2 = 20

    for r_id, (ref_acc, (seq, _)) in enumerate(help_functions.readfq(refs)):
        ref_id_to_accession[r_id] = ref_acc
        for positions, hash_val in seq_to_mixedaltstrobes_iter(seq, k_size1, k_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator):
            idx[hash_val].append(r_id)
            for pos in positions:
                idx[hash_val].append(pos)
            cntr += 1
            if cntr % 1000000 == 0:
                print("{0} mixedaltstrobes created from references".format(cntr))
        # print(hash_val, r_id, pos)
    return idx, ref_id_to_accession, cntr






def seq_to_altstrobes_iter(seq: str, k_size1: int, k_size2: int, strobe_w_min_offset: int,
                           strobe_w_max_offset: int, prime: int, w: int,
                           order: int) -> Iterator[tuple]:
    """
    Iterator for creation of altstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size1/2: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: an iterator for creating altstrobes
    """

    hash_seq_list1 = [(i, hash(seq[i:i+k_size1])) for i in range(len(seq) - k_size1 + 1)]
    hash_seq_list2 = [(i, hash(seq[i:i+k_size2])) for i in range(len(seq) - k_size2 + 1)]
    strobe_w_max_offset2 = strobe_w_max_offset - k_size1 # 40 vs 30
    strobe_w_min_offset2 = strobe_w_min_offset - k_size1 # 21 vs 11

    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned1 = thinner([h for i, h in hash_seq_list1], w)
        hash_seq_list_thinned2 = thinner([h for i, h in hash_seq_list2], w)
    else:
        hash_seq_list_thinned1 = hash_seq_list1
        hash_seq_list_thinned2 = hash_seq_list2

    for index_position, (p1, hash_m1) in enumerate(hash_seq_list_thinned1):  # [:-k_size]:
        if p1 >= len(hash_seq_list2) - (order-1)*k_size2:
            break
        # hash_m1 = hash_seq_list[p]

        if hash_m1 % 2 == 0: # first x, then 2x (e.g. 10-20)
            # print("10-20")
            if p1 + (order-1) * strobe_w_max_offset2 <= len(hash_seq_list2):
                windows = list()
                for window_order in range(1, order):
                    start = p1 + strobe_w_min_offset2 + (window_order-1) * strobe_w_max_offset2
                    end = min(p1 + window_order * strobe_w_max_offset2, len(hash_seq_list2))
                    windows.append((start, end))

            else:
                windows = list()
                for window_order in range(1, order):
                    start = (max(
                        p1+window_order*k_size2,
                        len(hash_seq_list2) + strobe_w_min_offset2 - (order - window_order) * strobe_w_max_offset2
                        )
                    )

                    end = min(p1 + window_order * strobe_w_max_offset2, len(hash_seq_list2))
                    windows.append((start, end))

            index = [p1, ] # -p1
            min_values = []
            min_hash_val = hash_m1
            for index_order in range(1, order):
                min_index, min_value = argmin([
                    (min_hash_val + hash_seq_list2[i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])

                min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list2[windows[index_order-1][0] + min_index][1]
                index.append(min_index+windows[index_order-1][0])
                index.append(min_index+windows[index_order-1][0]+k_size1)
                min_values.append(min_value)

        else: # first 2x, then x (e.g. 20-10)
            # print("20-10")
            if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list1):
                windows = list()
                for window_order in range(1, order):
                    start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                    end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list1))
                    windows.append((start, end))

            else:
                windows = list()
                for window_order in range(1, order):
                    start = (max(
                        p1+window_order*k_size1,
                        len(hash_seq_list1) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                        )
                    )

                    end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list1))
                    windows.append((start, end))

            index = [p1, p1+k_size1]
            min_values = []
            min_hash_val = hash_seq_list2[index_position][1]
            for index_order in range(1, order):
                min_index, min_value = argmin([
                    (min_hash_val + hash_seq_list1[i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])

                min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list1[windows[index_order-1][0] + min_index][1]
                index.append(min_index+windows[index_order-1][0]) # -
                min_values.append(min_value)

        yield index, min_hash_val


def altstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
               strobe_w_max_offset: int, w: int, order: int = 2,
               prime: int = 997) -> dict:
    """
    Strobemer seeding protocol to sample altstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :returns: a dictionary with positions along the string as keys and the altstrobes as value
    """

    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % (order+1) != 0:
        print("WARNING: kmer size {0} is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size1 = int(k_size/3)
    m_size2 = int(2*k_size/3)

    altstrobes = {tuple(index): h for index, h, min_values in seq_to_altstrobes_iter(
        seq, m_size1, m_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order
    )}
    return altstrobes

def altstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                     strobe_w_max_offset: int, w: int, order: int = 2,
                     buffer_size: int = 10000000) -> Iterator[tuple]:
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in altstrobes(
                substring, k_size, strobe_w_min_offset, strobe_w_max_offset,
                w, order=order).items():

            yield p, m

def seq_to_mixedaltstrobes_iter(seq: str, k_size1: int, k_size2: int, strobe_w_min_offset: int,
                                strobe_w_max_offset: int, prime: int, w: int,
                                order: int, denominator: int, numerator: int) -> Iterator[tuple]:
    """
    Iterator for creating of mixedaltstrobes of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: an iterator for creating mixedaltstrobes
    """
    hash_seq_list1 = [(i, hash(seq[i:i+k_size1])) for i in range(len(seq) - k_size1 + 1)]
    hash_seq_list2 = [(i, hash(seq[i:i+k_size2])) for i in range(len(seq) - k_size2 + 1)]
    strobe_w_max_offset2 = strobe_w_max_offset - k_size1 # 40 vs 30
    strobe_w_min_offset2 = strobe_w_min_offset - k_size1 # 21 vs 11

    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned1 = thinner([h for i, h in hash_seq_list1], w)
        hash_seq_list_thinned2 = thinner([h for i, h in hash_seq_list2], w)
    else:
        hash_seq_list_thinned1 = hash_seq_list1
        hash_seq_list_thinned2 = hash_seq_list2

    for index_position, (p1, hash_m1) in enumerate(hash_seq_list_thinned1):  # [:-k_size]:
        if p1 >= len(hash_seq_list2) - (order-1)*k_size2:
            break
        # hash_m1 = hash_seq_list[p]

        if hash_m1 % denominator < numerator:  # pick altstrobe
            if hash_m1 % 2 == 0: # first x, then 2x (e.g. 10-20)
                # print("10-20")
                if p1 + (order-1) * strobe_w_max_offset2 <= len(hash_seq_list2):
                    windows = list()
                    for window_order in range(1, order):
                        start = p1 + strobe_w_min_offset2 + (window_order-1) * strobe_w_max_offset2
                        end = min(p1 + window_order * strobe_w_max_offset2, len(hash_seq_list2))
                        windows.append((start, end))

                else:
                    windows = list()
                    for window_order in range(1, order):
                        start = (max(
                            p1+window_order*k_size2,
                            len(hash_seq_list2) + strobe_w_min_offset2 - (order - window_order) * strobe_w_max_offset2
                            )
                        )

                        end = min(p1 + window_order * strobe_w_max_offset2, len(hash_seq_list2))
                        windows.append((start, end))

                index = [p1, ] # -p1
                min_values = []
                min_hash_val = hash_m1
                for index_order in range(1, order):
                    min_index, min_value = argmin([
                        (min_hash_val + hash_seq_list2[i][1]) % prime
                        for i in range(*windows[index_order-1])
                    ])

                    min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list2[windows[index_order-1][0] + min_index][1]
                    index.append(min_index+windows[index_order-1][0])
                    index.append(min_index+windows[index_order-1][0]+k_size1)
                    min_values.append(min_value)

            else: # first 2x, then x (e.g. 20-10)
                # print("20-10")
                if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list1):
                    windows = list()
                    for window_order in range(1, order):
                        start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                        end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list1))
                        windows.append((start, end))

                else:
                    windows = list()
                    for window_order in range(1, order):
                        start = (max(
                            p1+window_order*k_size1,
                            len(hash_seq_list1) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                            )
                        )

                        end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list1))
                        windows.append((start, end))

                index = [p1, p1+k_size1]
                min_values = []
                min_hash_val = hash_seq_list2[index_position][1]
                for index_order in range(1, order):
                    min_index, min_value = argmin([
                        (min_hash_val + hash_seq_list1[i][1]) % prime
                        for i in range(*windows[index_order-1])
                    ])

                    min_hash_val = min_hash_val + (index_order * (-1)**index_order) * hash_seq_list1[windows[index_order-1][0] + min_index][1]
                    index.append(min_index+windows[index_order-1][0]) # -

            yield index, min_hash_val

        else:  # pick k-mer
            index = tuple(p1 + (strobe_num) * k_size1 for strobe_num in range(order+1))
            l = int(k_size1+k_size2/2*order)
            yield index, hash(seq[p1: p1+l])


def mixedaltstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int, order: int = 2,
                    strobe_fraction: float = 0.5) -> dict:
    """
    Strobemer seeding protocol to sample mixedaltstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :returns: a dictionary with positions along the string as keys and the altstrobes as value
    """
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % (order+1) != 0:
        print("WARNING: kmer size is not evenly divisible with {0}, will use {1} as kmer size: ".format(order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size1 = int(k_size/3)
    m_size2 = int(2*k_size/3)

    mixedaltstrobes = {tuple(index): h for index, h in seq_to_mixedaltstrobes_iter(
        seq, m_size1, m_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator
    )}

    return mixedaltstrobes

# def sort_merge(sorted_list):
#     sort_merged_list = []
#     curr_merge = sorted_list[0]
#     for i, t1 in enumerate( sorted_list[:-1] ):
#         r_id, r_pos, q_pos, length = t1
#         r2_id, r2_pos, q2_pos, length2 = sorted_list[i+1]
#         # print(i, r_id, r_pos, r2_id, r2_pos)
#         # print(r2_pos, q2_pos)
#         if r_id == r2_id:
#             # print("OK", q2_pos <= q_pos + length <= q2_pos+ length, r2_pos <= r_pos + length <= r2_pos + length)
#             # print("2", q2_pos, q_pos + length, q2_pos+ length, r2_pos, r_pos + length, r2_pos + length)
#             # overlapping on both query and ref
#             # print(q2_pos + length2, q_pos + length, curr_merge[3])
#             if q2_pos <= q_pos + length <= q2_pos+ length  and r2_pos <= r_pos + length <= r2_pos + length:
#                 # curr_merge = (r_id, curr_merge[1], curr_merge[2], max(q2_pos + length2, q_pos + length ) -  q_pos ) # hit length on query sequence
#                 curr_merge = (r_id, curr_merge[1], curr_merge[2], max(r2_pos + length2, r_pos + length ) -  r_pos ) # hit length on reference sequence
#                 # print("HERER")

#             else:
#                 # time to add old element
#                 sort_merged_list.append(curr_merge)
#                 curr_merge = sorted_list[i+1]

#         else:
#             # time to add old element
#             sort_merged_list.append(curr_merge)
#             curr_merge = sorted_list[i+1]
#         # print(curr_merge)
#         # print(sort_merged_list)
#     # print(curr_merge)
#     sort_merged_list.append(curr_merge)
#     return sort_merged_list


def matches_query(strobes:list, idx: dict, k: int) -> Iterator:
    """
    Iterate over query in ascending order and check if hash value in reference

    :param strobes: list of kmers/strobemers with (positions, hash_values) for each sampled kmer/strobemer
    :param idx: a dictionary with reference indexes
    :param k: strobe size
    :returns: iterator with tuple of q_start, q_end and hash value
    """
    for q_positions, h in strobes:
        if h in idx:
            # extract query start and end position
            q_start = q_positions[0] + 1
            q_end = q_positions[-1] + k + 1
            yield q_start, q_end, h


def extend_nams(cpm: dict, merged_matches: list, r_id: int, q_start: int, q_end: int, r_start: int, r_end: int):
    """
    Extend Non-Overlapping Approximate Matches and updates the cmp dictionary.

    :param cpm: dict with open merge intervals
    :param merged_matches: list of merged matches
    :param r_id: reference id
    :param q_start: query start position of the current hit
    :param q_end: query end position of the current hit
    :param r_start: reference start position matching the query hit position
    :param r_end: reference end position matching the query hit position
    :returns: updated dictionary of open merge intervals and list of merged_matches
    """
    # query and reference positions are already in cmp
    # check all queries for overlap
    is_added_to_an_interval_query = False
    for prev_q_end in list(cpm[r_id].keys()):

        if q_start <= prev_q_end:  # overlap on query
            is_added_to_an_interval_query = True
            is_added_to_an_interval_reference = False

            # check all references for overlap
            for prev_r_end in list(cpm[r_id][prev_q_end].keys()):
                prev_q_start, prev_q_end, prev_r_start, prev_r_end = cpm[r_id][prev_q_end][prev_r_end]

                # update query positions
                new_q_end = max(prev_q_end, q_end)

                if prev_r_start <= r_start <= prev_r_end:  # overlap on reference
                    is_added_to_an_interval_reference = True
                    # update reference positions
                    new_r_end = max(prev_r_end, r_end)

                    # update cpm dictionary
                    del cpm[r_id][prev_q_end][prev_r_end]  # delete old reference positions
                    if not cpm[r_id][prev_q_end]:  # delete old query positions if empty
                        del cpm[r_id][prev_q_end]
                    if new_q_end not in cpm[r_id]:  # add new query and reference positions
                        cpm[r_id][new_q_end] = {}
                        cpm[r_id][new_q_end][new_r_end] = (prev_q_start, new_q_end, prev_r_start, new_r_end)
                    elif new_r_end not in cpm[r_id][new_q_end]:  # append
                        cpm[r_id][new_q_end][new_r_end] = (prev_q_start, new_q_end, prev_r_start, new_r_end)
                    else:  # was already present
                        (old_q_start, new_q_end, old_r_start, new_r_end) = cpm[r_id][new_q_end][new_r_end]
                        cpm[r_id][new_q_end][new_r_end] = (min(old_q_start, prev_q_start), new_q_end, min(old_r_start, prev_r_start), new_r_end)

                # no overlap with any reference
                if not is_added_to_an_interval_reference:
                    if new_q_end not in cpm[r_id]:  # new added 1
                        cpm[r_id][new_q_end] = {}
                        cpm[r_id][new_q_end][r_end] = (q_start, new_q_end, r_start, r_end)

                    elif r_end not in cpm[r_id][new_q_end]:  # new added 2
                        cpm[r_id][new_q_end][r_end] = (q_start, new_q_end, r_start, r_end)

                    else:  # was already present
                        (old_q_start, new_q_end, old_r_start, new_r_end) = cpm[r_id][new_q_end][r_end]
                        cpm[r_id][new_q_end][new_r_end] = (min(old_q_start, q_start), new_q_end, min(old_r_start, r_start), new_r_end)

        else:
            # revove the intervals that we have passed on the query here to not make the cpm dict too large...
            # add to merged_matches dict
            for r_pos_stop in cpm[r_id][prev_q_end]:
                (q_pos_start, q_pos_stop, r_pos, r_pos_stop) = cpm[r_id][prev_q_end][r_pos_stop]
                merged_matches.append((r_id, r_pos, q_pos_start, r_pos_stop - r_pos))
            del cpm[r_id][prev_q_end]
            # add match to cpm if all all intervals have been passed
            if q_end not in cpm[r_id]:
                cpm[r_id][q_end] = {r_end: (q_start, q_end, r_start, r_end)}

    if not is_added_to_an_interval_query:  # no overlap with prev query sequences
        cpm[r_id][q_end] = {}
        cpm[r_id][q_end][r_end] = (q_start, q_end, r_start, r_end)

    return cpm, merged_matches


def get_unmerged_matches(strobes: list, idx: dict, k: int,
                         ref_id_to_accession: dict, acc, selfalign: bool,
                         order: int) -> list:
    """
    It is seriously advised to merge matches as the files can become huge otherwise and fill up all diskspace.

    :param strobes: list of kmers/strobemers with (positions, hash_values) for each sampled kmer/strobemer
    :param idx: a dictionary with reference indexes
    :param k: strobe size
    :param ref_id_to_accession: dictionary with references and accession numbers
    :param acc: accession number
    :param selfalign: aligns sequences to itself (mainly for bugfixing)
    :param order: number of substrings/strobes
    :returns: sorted list of (merged) matches
    """
    matches = []
    # iterate over query in ascending order and check if hash value in reference
    for q_start, q_end, h in matches_query(strobes, idx, k):
        for r_id, *r_positions in grouper(idx[h], order+1):
            matches.append((r_id, r_positions[0] + 1, q_start, r_positions[-1] - r_positions[0] + k))
    return sorted(matches, key=lambda x: (x[0], x[2], x[1]))


def get_unmerged_matches_altstrobes(strobes: list, idx: dict, k: int,
                         ref_id_to_accession: dict, acc, selfalign: bool,
                         order: int) -> list:
    """
    It is seriously advised to merge matches as the files can become huge otherwise and fill up all diskspace.

    :param strobes: list of kmers/strobemers with (positions, hash_values) for each sampled kmer/strobemer
    :param idx: a dictionary with reference indexes
    :param k: strobe size
    :param ref_id_to_accession: dictionary with references and accession numbers
    :param acc: accession number
    :param selfalign: aligns sequences to itself (mainly for bugfixing)
    :param order: number of substrings/strobes
    :returns: sorted list of (merged) matches
    """
    matches = []
    # iterate over query in ascending order and check if hash value in reference
    for q_start, q_end, h in matches_query(strobes, idx, k):
        for r_id, *r_positions in grouper(idx[h], order+1):
            # matches.append((r_id, r_positions[0] + 1, q_start, r_positions[-1] - r_positions[0] + k))
            matches.append((r_id, r_positions[0] + 1, q_start, r_positions[1] - r_positions[0] + k))
            matches.append((r_id, r_positions[1] + 1, q_start + r_positions[1]-r_positions[0], r_positions[-1] - r_positions[1] + k))
    return sorted(matches, key=lambda x: (x[0], x[2], x[1]))


def get_merged_matches(strobes: list, idx: dict, k: int, ref_id_to_accession: dict,
                       acc, selfalign: bool, order: int) -> list:
    """
    The merging of matches is a simple linear merging. If there are repetitive matches across e.g. a chromosome
    the merging will be broken up at the repetitive kmer. To solve the merging exactly, we would need
    to solve the collinear chaining problem after we have out matches. There is no such functionality here.

    Another way to solve this is to do a post merging after sorting the merged matches.
    If two merged matches also overlaps, they can be merged again.

    :param strobes: list of kmers/strobemers with (positions, hash_values) for each sampled kmer/strobemer
    :param idx: a dictionary with reference indexes
    :param k: strobe size
    :param ref_id_to_accession: dictionary with references and accession numbers
    :param acc: accession number
    :param selfalign: aligns sequences to itself (mainly for bugfixing)
    :param order: number of substrings/strobes
    :returns: sorted list of (merged) matches
    """
    cpm = {}  # current potential merges
    # FORMAT: cpm[r_id] = {end_q: {end_r: (start_q, end_q, start_r, end_r)}}
    merged_matches = []

    # iterate over query in ascending order and check if hash value in reference
    for q_start, q_end, h in matches_query(strobes, idx, k):

        # iterate over references, all in ascending order
        for r_id, *r_positions in grouper(idx[h], order+1):

            # remove self matches with below if statement
            if not selfalign and ref_id_to_accession[r_id] == acc:
                continue

            # extract reference start and end position
            r_start = r_positions[0] + 1
            r_end = r_positions[-1] + k + 1

            if r_id in cpm:  # extend non-overlapping approximate matches and update merged_matches and cpm
                cpm, merged_matches = extend_nams(cpm, merged_matches, r_id, q_start, q_end, r_start, r_end)

            else:
                # add query and reference positions (start+end) to cpm
                cpm[r_id] = {q_end: {r_end: (q_start, q_end, r_start, r_end)}}
                continue

    # close all open merge intervals
    for r_id in cpm.keys():
        for q_stop in cpm[r_id]:
            for r_stop in cpm[r_id][q_stop]:
                (q_p1, q_pos_stop, r_pos, r_pos_stop) = cpm[r_id][q_stop][r_stop]
                merged_matches.append((r_id, r_pos, q_p1, r_pos_stop - r_pos))

    # if no matches were found, an empty list is returned
    if not merged_matches:
        return []

    # remove duplicated matches
    unique_merged_matches = dict()
    for match in merged_matches:
        r_id, r_pos, q_pos, r_len = match
        # add match to unique matches if no match with same query and reference ending
        if not (r_pos+r_len, q_pos+r_len) in unique_merged_matches:
            unique_merged_matches[(r_pos+r_len, q_pos+r_len)] = (r_id, r_pos, q_pos, r_len)
        # replace unique match by longer match with same query and reference ending
        elif r_len > unique_merged_matches[(r_pos+r_len, q_pos+r_len)][3]:
            unique_merged_matches[(r_pos+r_len, q_pos+r_len)] = (r_id, r_pos, q_pos, r_len)
        # discard current match as a longer one is already in unique_merged_matches
        else:
            continue

    # return sorted(set(merged_matches), key=lambda x: (x[0], x[2], x[1]))
    return sorted(set(unique_merged_matches.values()), key=lambda x: (x[0], x[3], x[2], x[1]), reverse=True)


def print_matches_to_file(query_matches: list, ref_id_to_accession: dict,
                          outfile: str, reverse: bool) -> None:
    """
    Print matches to tsv-file.

    Output format is a tab separated file on the same format as MUMmer,
    with identical fields except the last one which is approximate reference
    sequence match length instead of what MUMmer produce:

    >query_accession
        ref_id  ref_pos query_pos   match_length_on_reference

    :param query_matches: list of (merged) query matches
    :param ref_id_to_accession: dictionary with references and accession numbers
    :param outfile: string with outfile path and name
    :param reverse: a bool indicating whether reverse complement were matched
    """
    for q_acc, read_matches in query_matches:
        if reverse:
            outfile.write("> {0} Reverse\n".format(q_acc))
        else:
            outfile.write("> {0}\n".format(q_acc))
        for (r_id, ref_p, q_pos, k) in read_matches:
                ref_acc = ref_id_to_accession[r_id]
                outfile.write("  {0} {1} {2} {3}\n".format(ref_acc, ref_p, q_pos, k))


def main(args):
    """
    """
    PRIME = 997
    w = args.w

    if not args.kmer_index:
        print("\nUsing")
        print(f'Order: {args.n}')
        print(f'k: {args.k}')
        print(f'w_min: {args.strobe_w_min_offset}')
        print(f'w_max: {args.strobe_w_max_offset}')
        print("\n")
    else:
        print("\nUsing")
        print("kmers")
        print(f'k: {args.k}')



    if args.strobe_fraction:
        fraction = Fraction(str(args.strobe_fraction))
        denominator = fraction.denominator
        numerator = fraction.numerator

    if args.kmer_index:
        args.n = 1
        idx, ref_id_to_accession, cntr = build_kmer_index(open(args.references, 'r'), args.k, w)
        print("{0} kmers created from references\n".format(cntr))
        # print(idx)
    elif args.minstrobe_index:
        idx, ref_id_to_accession, cntr = build_strobemer_index("minstrobes", open(args.references, 'r'), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)
        print("{0} minstrobes created from references\n".format(cntr))
    elif args.randstrobe_index:
        idx, ref_id_to_accession, cntr = build_strobemer_index("randstrobes", open(args.references, 'r'), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)
        print("{0} randstrobes created from references\n".format(cntr))
    elif args.hybridstrobe_index:
        idx, ref_id_to_accession, cntr = build_hybridstrobe_index(open(args.references, 'r'), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)
        print("{0} hybridstrobes created from references\n".format(cntr))
    elif args.altstrobe_index:
        args.k = 10
        idx, ref_id_to_accession, cntr = build_altstrobe_index(open(args.references, 'r'), 10, 20, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)
        print("{0} altstrobes created from references\n".format(cntr))
    elif args.mixedminstrobe_index:
        idx, ref_id_to_accession, cntr = build_mixedstrobemer_index("mixedminstrobes", open(args.references, 'r'), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)
        print("{0} mixedminstrobes created from references\n".format(cntr))
    elif args.mixedrandstrobe_index:
        idx, ref_id_to_accession, cntr = build_mixedstrobemer_index("mixedrandstrobes", open(args.references, 'r'), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)
        print("{0} mixedrandstrobes created from references\n".format(cntr))
    elif args.mixedhybridstrobe_index:
        idx, ref_id_to_accession, cntr = build_mixedhybridstrobe_index(open(args.references, 'r'), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)
        print("{0} mixedhybridstrobes created from references\n".format(cntr))
    elif args.mixedaltstrobe_index:
        args.k = 10
        idx, ref_id_to_accession, cntr = build_mixedaltstrobe_index(open(args.references, 'r'), 10, 20, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)
        print("{0} mixedaltstrobes created from references\n".format(cntr))
    else:
        raise NameError

    outfile = open(os.path.join(args.outfolder, args.prefix + ".txt"), 'w')
    outfile_summary = open(os.path.join(args.outfolder, args.prefix + " (Summary).txt"), 'w')
    print("\nProcessing query sequences")
    for n_query_sequences, (acc, (full_seq, _)) in enumerate(help_functions.readfq(open(args.queries, 'r'))):
        print("Processing query sequence #", n_query_sequences, "of length", len(full_seq))
        results = {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0, "coll": 0}

        for segment in range(math.ceil(len(full_seq)/2000)):
            seq = full_seq[segment*2000:segment*2000+2000]
            tmp_outfile = open("tmp.tsv", 'w')
            if len(seq) > args.k + args.strobe_w_max_offset:

                query_matches = []

                if args.rev_comp:
                    matches_rc = []

                if args.kmer_index:
                    strobes = [(positions, h) for positions, h in seq_to_kmer_iter(seq, args.k, w)]
                    read_matches = get_unmerged_matches(strobes, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                elif args.minstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_minstrobes_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                    read_matches = get_unmerged_matches(strobes, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                elif args.randstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_randstrobes_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                    read_matches = get_unmerged_matches(strobes, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                elif args.hybridstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_hybridstrobes_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                    read_matches = get_unmerged_matches(strobes, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                elif args.altstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_altstrobes_iter(seq, 10, 20, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                    read_matches = get_unmerged_matches_altstrobes(strobes, idx, 10, ref_id_to_accession, acc, args.selfalign, args.n+1)

                elif args.mixedminstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_mixedminstrobes_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                    read_matches = get_unmerged_matches(strobes, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                elif args.mixedrandstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_mixedrandstrobes_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                    read_matches = get_unmerged_matches(strobes, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                elif args.mixedhybridstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_mixedhybridstrobes_iter(seq, args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                    read_matches = get_unmerged_matches(strobes, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                elif args.mixedaltstrobe_index:
                    strobes = [(positions, h) for positions, h in seq_to_mixedaltstrobes_iter(seq, 10, 20, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                    read_matches = get_unmerged_matches_altstrobes(strobes, idx, 10, ref_id_to_accession, acc, args.selfalign, args.n+1)

                else:
                    raise NameError

                # print(strobes)
                query_matches.append((acc, read_matches))
                print_matches_to_file(query_matches, ref_id_to_accession, tmp_outfile, False)

                if args.rev_comp:
                    if args.kmer_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_kmer_iter(reverse_complement(seq), args.k, w)]
                        read_matches_rc = get_unmerged_matches(strobes_rc, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                    elif args.minstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_minstrobes_iter(reverse_complement(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                        read_matches_rc = get_unmerged_matches(strobes_rc, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                    elif args.randstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_randstrobes_iter(reverse_complement(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                        read_matches_rc = get_unmerged_matches(strobes_rc, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                    elif args.hybridstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_hybridstrobes_iter(reverse_complement(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                        read_matches_rc = get_unmerged_matches(strobes_rc, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                    elif args.altstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_altstrobes_iter(reverse_complement(seq), 10, 20, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n)]
                        read_matches_rc = get_unmerged_matches_altstrobes(strobes_rc, idx, 10, ref_id_to_accession, acc, args.selfalign, args.n+1)

                    elif args.mixedminstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_mixedminstrobes_iter(reverse_complement(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                        read_matches_rc = get_unmerged_matches(strobes_rc, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                    elif args.mixedrandstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_mixedrandstrobes_iter(reverse_complement(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                        read_matches_rc = get_unmerged_matches(strobes_rc, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                    elif args.mixedhybridstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_mixedhybridstrobes_iter(reverse_complement(seq), args.k, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                        read_matches_rc = get_unmerged_matches(strobes_rc, idx, args.k, ref_id_to_accession, acc, args.selfalign, args.n)

                    elif args.mixedaltstrobe_index:
                        strobes_rc = [(positions, h) for positions, h in seq_to_mixedaltstrobes_iter(reverse_complement(seq), 10, 20, args.strobe_w_min_offset, args.strobe_w_max_offset, PRIME, w, args.n, denominator, numerator)]
                        read_matches_rc = get_unmerged_matches_altstrobes(strobes_rc, idx, 10, ref_id_to_accession, acc, args.selfalign, args.n+1)

                    else:
                        raise NameError

                    matches_rc.append((acc, read_matches_rc))
                    print_matches_to_file(matches_rc, ref_id_to_accession, tmp_outfile, True)
                tmp_outfile.close()

                read_coverage_solution = {}
                total_disjoint_matches = 0
                tot_genome_length = 0

                for (query_acc, nams) in get_NAM_records("tmp.tsv", acc):
                    q_acc = query_acc.split()[0]
                    tot_genome_length += len(seq)
                    # print(t_chr_id, t_start, t_end)
                    for ref_id in nams:
                        chrom_nams = nams[ref_id]
                        total_disjoint_matches += len(chrom_nams)
                        solutions, opt_cov = n_logn_read_coverage(chrom_nams)
                        # pick best from forward and reverse strand
                        if q_acc in read_coverage_solution:
                            c_prev, _ = read_coverage_solution[q_acc]
                            if c_prev < opt_cov:
                                # print(ref_id, opt_cov)
                                read_coverage_solution[q_acc] = (opt_cov, solutions[0])
                        else:
                            # print(ref_id, opt_cov)
                            read_coverage_solution[q_acc] = (opt_cov, solutions[0])

                tot_genome_length = tot_genome_length/2 # remove double counting of reverse complements
                collinear_chain_nam_sizes = []
                total_bp_covered = 0

                sc_positions = []
                gaps = []
                gap_pos = 0

                # collinear_outfile = open("collinear_matches_out.tsv", "w")
                for q_acc in read_coverage_solution:
                    # collinear_outfile.write("> {0}\n".format(q_acc))
                    opt_cov, solution = read_coverage_solution[q_acc]
                    total_bp_covered += opt_cov
                    for n in solution:
                        collinear_chain_nam_sizes.append(n.y - n.x)
                        sc_positions.append(n.c)
                        sc_positions.append(n.c + n.val - args.k)

                        if n.c > gap_pos:
                            gaps.append(n.c-gap_pos-1)
                        gap_pos = n.d
                        # collinear_outfile.write("  {0} {1} {2} {3}\n".format(n.chr_id, n.x, n.c, n.val))

                gaps.append(len(seq)-gap_pos)
                coll_esize = e_size(collinear_chain_nam_sizes, tot_genome_length)
                # collinear_outfile.close()

                m = len(n)
                #m = len(read_matches) + len(read_matches_rc)
                if args.altstrobe_index or args.mixedaltstrobe_index:
                    m = m/2 # correct for representation as two hits
                # sc = get_sequence_coverage(sorted(sc_positions), args.k)
                sc = get_sequence_coverage(sorted(sc_positions), args.k)
                seeds = len(strobes) + len(strobes_rc)

                if seeds > 1:
                    outfile.write("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9}\n".format(
                        n_query_sequences+1, # read number
                        len(full_seq), # read length
                        segment+1, # segment number
                        len(seq), # sequence length
                        total_disjoint_matches, # m
                        round(m/seeds*2, 2), # mp
                        round(sc/len(seq), 2), # sc
                        round(total_bp_covered/tot_genome_length,2), # mc
                        round(get_e_size(gaps, len(seq), 1),2), # gaps
                        round(coll_esize/2,1) #coll_esize
                    ))

                    results["m"] += m
                    results["mp"] += seeds
                    results["sc"] += sc
                    results["mc"] += total_bp_covered
                    results["gaps"].append(gaps)
                    results["coll"] += coll_esize/2


        flat = [g for l in results["gaps"] for g in l]
        if flat:
            # avg_island_len = sum(flat)/len(flat)
            # print(protocol)
            q_e_size = get_e_size(flat, len(full_seq), 1)
        # else:
        #     avg_island_len = 0
        res = [
            round(100*results["m"]*2/results["mp"], 1), # account for the fact that a seed can only match one direction (see collinear solution)
            100*results["sc"]/len(full_seq),
            100*results["mc"]/len(full_seq),
            q_e_size,
            results["coll"]
        ]
        outfile_summary.write(
        ">Query " + str(n_query_sequences) + " & " + str(len(full_seq)) + " & " + " & ".join([str(round(r, 1)) for r in res])+"\n")
    outfile.close()
    outfile_summary.close()

        # print("Finished processing {0} query sequences.".format(n_query_sequences+1))

def get_sequence_coverage(positions, k_len):
    covered_bases = 0

    if len(positions) == 0:
        return 0

    prev_p = positions[0]
    covered_bases += k_len  # for first pos
    if len(positions) == 1:
        return covered_bases

    for p in positions[1:]:
        if p <= prev_p + k_len - 1:
            covered_bases += p-prev_p
        else:
            covered_bases += k_len
        prev_p = p
    return covered_bases


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--queries', type=str,  default=False, help='Path to query fasta or fastq file')
    parser.add_argument('--references', type=str,  default=False, help='Path to reference fasta or fastq file')
    parser.add_argument('--k', type=int, default=15, help='Length of kmer/all strobes combined')
    parser.add_argument('--strobe_w_min_offset', type=int, default=20, help='Strobemer window start offset from first k-mer. If kmer start at pos i, first\
                                                                            window will start at i+strobe_w_min_offset. Default: 20nt donwstream from start of first kmer.')
    parser.add_argument('--strobe_w_max_offset', type=int, default=70, help='Strobemer window end. If kmer start at pos i, first\
                                                                            window will stop at i+strobe_w_max_offset. Default: 70nt donwstream from start of first kmer.')
    parser.add_argument('--w', type=int, default=1, help='Thinning window size applied to reference sequences (default = 1, i.e., no thinning)')
    parser.add_argument('--n', type=int, default=2, help='Order on strobes')
    parser.add_argument('--strobe_fraction', type=float, default=1, help='Fraction of sampled strobemers, rest kmers')
    parser.add_argument('--dont_merge_matches', action="store_true",  help='Do not merge matches with this option. It is seriously advised to\
                                                                     merge matches as the files can become huge otherwise and fill up all diskspace.\
                                                                     Do not specify this option unless you know what you are doing! Mostly here for\
                                                                     development/bugchecking purposas. The default option is to merge matches if they\
                                                                     are consectutive on both query and reference to create MAM-like matches \
                                                                     (maximal approximate matches) of various lengths, much like the output of MUMmer. This is\
                                                                     disk space frendilier, although these files can get large too.')
    parser.add_argument('--outfolder', type=str,  default=None, help='Folder to output TSV match file.')
    parser.add_argument('--prefix', type=str,  default="matches", help='Filename prefix (default "matches").')
    parser.add_argument('--kmer_index', action="store_true",  help='Produce non-overlapping approximate matches for k-mers')
    parser.add_argument('--minstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for minstrobes')
    parser.add_argument('--randstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for randstrobes')
    parser.add_argument('--hybridstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for hybridstrobes')
    parser.add_argument('--altstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for altstrobes')
    parser.add_argument('--mixedminstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for mixed minstrobes/kmers based on --strobe_fraction')
    parser.add_argument('--mixedrandstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for mixed randstrobes/kmers based on --strobe_fraction')
    parser.add_argument('--mixedhybridstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for mixed hybridstrobes/kmers based on --strobe_fraction')
    parser.add_argument('--mixedaltstrobe_index', action="store_true",  help='Produce non-overlapping approximate matches for mixed altstrobes/kmers based on --strobe_fraction')
    parser.add_argument('--selfalign', action="store_true",  help='Aligns sequences to itself (mainly used for bugfixing). Default is not align\
                                                                    sequences to themselves if the same file is given as references and queries.')
    # parser.add_argument('--compress', type=str,  default=None, help='Compress output')
    parser.add_argument('--rev_comp', action="store_true",  help='Match reverse complement of reads (output to separate file)')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    # if args.w != 1:
    #     raise NotImplementedError("Currently only w=1 is allowed, i.e., no thinning is implemented")

    main(args)