#!/usr/bin/env python3.9
# -*- coding: utf-8 -*

"""Short Description.

Long Description
"""

__authors__ = ["Benjamin D. Maier"]
__copyright__ = "Copyright Benjamin D. Maier & Kristoffer Sahlin | Sahlin Group"
__organization__ = "Department of Mathematics, Science for Life Laboratory, Stockholm University, 106 91, Stockholm, Sweden."
__credits__ = ["Benjamin D. Maier & Kristoffer Sahlin"]
__contact__ = "bmaier [at] ebi.ac.uk"
__date__ = "2023/03/10"
__created__ = "2022/02/XX"
__deprecated__ = False
__license__ = "MIT"
__maintainer__ = "Kristoffer Sahlin"
__email__ = "kristoffer.sahlin [at] scilifelab.se"
__status__ = "DSML Lvl. 1 - Concept"

import sys
import operator
import random
import copy
from collections import defaultdict, deque
from fractions import Fraction
from typing import Iterator
from itertools import chain

BITS = sys.hash_info.width
MAX = sys.maxsize
MAX_HASH_VALUE = int((2**BITS)/2) - 1


def Sequence(L: int) -> str:
    """
    Generate a random canonical DNA sequence.

    :param L: an integer representing the desired sequence length
    :returns: a string with a random DNA sequence of length L
    """
    DNA = "".join([random.choice("ACGT") for i in range(L)])
    return DNA


def argmin(array: list) -> tuple:
    """
    Find the value of x which minimizes f(x) over the set of candidates for x

    :param array: a list to minimize
    :returns: a tuple with the index position and the value of the lowest element
    """
    min_index = array.index(min(array))
    min_val = array[min_index]
    return min_index, min_val


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


def minimizers(seq: str, k_size: int, w_size: int) -> list:
    """
    Sample the smallest k-mer by hash-value in a pre-defined ordering of each k-mer in the window

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the minimizers
    :param w_size: window size
    :returns: a list with minimizers
    """
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w+1)])
    curr_min = min(window_kmers)
    minimizers = [(curr_min, list(window_kmers).index(curr_min))]

    for i in range(w+1, len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer:
            curr_min = min(window_kmers)
            minimizers.append((curr_min, list(window_kmers).index(curr_min) + i - w))

        # previous minimizer still in window, we only need to compare with the recently added kmer
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append((curr_min, i))

    return minimizers


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


def kmers(seq: str, k_size: int, w: int) -> dict:
    """
    Sample a substrings of length k contained within a given sequence

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: a dictionary with positions along the string as keys and the kmers as value
    """
    if w > 1:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
        kmers_pos = {i: h for i, h in hash_seq_list_thinned}
    else:
        kmers_pos = {i: hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)}

    return kmers_pos


def kmer_iter(seq: str, k_size: int, w: int) -> Iterator[tuple]:
    """
    Iterator for creating kmers

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: an iterator for creating kmers
    """
    if w > 1:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
        for p, h in hash_seq_list_thinned:
            yield p, h
    else:
        hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
        for p, h in hash_seq_list:
            yield p, h


def spaced_kmers(seq: str, k_size: int, span_size: int, positions: set,
                 w: int) -> dict:
    """
    Sample kmers with spaced seeds

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the strobe
    :param span_size: length between first and last position
    :param positions: a set of positions to consider for the spaced k-mer
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :returns: a dictionary with positions along the string as keys and the spaced kmer as value
    """
    assert len(positions) == k_size

    if w > 1:
        hash_seq_list = [
            (i, hash("".join([seq[i + j] for j in range(span_size) if j in positions])))
            for i in range(len(seq) - span_size + 1)
        ]
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
        spaced_kmers_pos = {i: h for i, h in hash_seq_list_thinned}
    else:
        spaced_kmers_pos = {
            i: hash("".join([seq[i + j] for j in range(span_size) if j in positions]))
            for i in range(len(seq) - span_size + 1)
        }
    # print(positions, len(positions), span_size)
    # well, this is not the most time efficient way to sample spaced kmers but works for simulations...
    return spaced_kmers_pos


def spaced_kmers_iter(seq: str, k_size: int, span_size: int,
                      positions: set) -> Iterator[int]:
    """
    Iterator for generating spaced kmers

    :param seq: a string with a nucleotide sequence
    :param k_size: length of the kmers
    :param span_size: length between first and last position
    :param positions: a set of positions to consider for the spaced k-mer
    :returns: an iterator for creating spaced_kmers
    """
    assert len(positions) == k_size
    # print(positions, len(positions), span_size)
    # well, this is not the most time efficient way to sample spaced kmers but works for simulations...
    for i in range(len(seq) - span_size + 1):
        yield hash("".join([seq[i + j] for j in range(span_size) if j in positions]))


def randstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                strobe_w_max_offset: int, w: int, order: int = 2,
                prime: int = 997, return_min_value: bool = False) -> dict:
    """
    Strobemer seeding protocol to sample randstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param return_min_value: a bool to specify whether the hash value we minimize in the function deciding how to pick the strobe should be returned
    :returns: a dictionary with positions along the string as keys and the randstrobes as value
    """
    randstrobes = dict()
    randstrobes_hash = dict()
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order

    if return_min_value:
        randstrobes = dict()
        randstrobes_min_values = dict()
        for index, h in seq_to_randstrobes_iter(seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order):
            randstrobes[tuple(index)] = h
            randstrobes_min_values[h] = min_values
        return randstrobes, randstrobes_min_values

    else:
        randstrobes = {tuple(index): h for index, h in seq_to_randstrobes_iter(
            seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order
        )}
        return randstrobes


def randstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                     strobe_w_max_offset: int, w: int, order: int = 2,
                     buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating randstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating randstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in randstrobes(
                substring, k_size, strobe_w_min_offset, strobe_w_max_offset,
                w, order=order).items():

            yield p, m


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
    prime = 997
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 >= len(hash_seq_list) - (order-1)*k_size:
            break
        # hash_m1 = hash_seq_list[p]

        if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list):
            windows = list()
            for window_order in range(1, order):
                start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                windows.append((start, end))

        else:
            windows = list()
            for window_order in range(1, order):
                start = (max(
                    p1+window_order*k_size,
                    len(hash_seq_list) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                    )
                )

                end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                windows.append((start, end))

        index = [p1, ]
        min_values = []
        min_hash_val = hash_m1
        for index_order in range(1, order):
            min_index, min_value = argmin([
                (hash_m1 ^ hash_seq_list[i][1])
                for i in range(*windows[index_order-1])
            ])

            min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list[windows[index_order-1][0] + min_index][1]
            index.append(min_index+windows[index_order-1][0])
            min_values.append(min_value)

        yield index, min_hash_val #, min_values


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
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
    # thinning
    if w > 1:
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 >= len(hash_seq_list) - (order-1)*k_size:
            break
        # hash_m1 = hash_seq_list[p]

        if hash_m1 % denominator < numerator:  # pick randstrobe
            if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list):
                windows = list()
                for window_order in range(1, order):
                    start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
                    end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                    windows.append((start, end))

            else:
                windows = list()
                for window_order in range(1, order):
                    start = (max(
                        p1+window_order*k_size,
                        len(hash_seq_list) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                        )
                    )

                    end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
                    windows.append((start, end))

            index = [p1, ]
            min_hash_val = hash_m1
            for index_order in range(1, order):
                min_index, min_value = argmin([
                    (hash_m1 ^ hash_seq_list[i][1])
                    for i in range(*windows[index_order-1])
                ])

                min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list[windows[index_order-1][0] + min_index][1]
                index.append(min_index+windows[index_order-1][0])
            yield index, min_hash_val

        else:  # pick k-mer
            index = tuple(p1 + (strobe_num) * k_size for strobe_num in range(order))
            yield index, hash(seq[p1: p1+k_size*order])


def mixedrandstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                     strobe_w_max_offset: int, w: int,
                     order: int = 2, strobe_fraction: float = 0.5) -> dict:
    """
    Mixed protocol to produce specified randstrobes and k-mers content fractions

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: a dictionary with positions along the string as keys and the mixed randstrobes/kmers as value
    """
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order

    mixedrandstrobes = {
        tuple(index): h
        for index, h in seq_to_mixedrandstrobes_iter(
            seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator
        )
    }

    return mixedrandstrobes


def mixedrandstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, w: int, order: int = 2,
                         strobe_fraction: float = 1,
                         buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedrandstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating mixedrandstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in mixedrandstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order, strobe_fraction).items():
            yield p, m


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
    hash_seq_list = [(i, hash(seq[i: i+k_size])) for i in range(len(seq) - k_size + 1)]

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


def minstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
               strobe_w_max_offset: int, w: int, order: int = 2) -> dict:
    """
    Strobemer seeding protocol to sample minstrobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: a dictionary with positions along the string as keys and the minstrobes as value
    """
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order
    assert m_size + (order-1) * strobe_w_max_offset < len(seq), "Last minstrobes window position is exceeding the sequence length, consider to use a lower order"

    minstrobes = {
        tuple(positions): h
        for positions, h in seq_to_minstrobes_iter(
            seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime,
            w, order
        )
    }

    return minstrobes


def minstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int, order: int = 2,
                    buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating minstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating minstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in minstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order=order).items():
            yield p, m


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
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
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
            index = tuple(p1 + (strobe_num) * k_size for strobe_num in range(order))
            yield index, hash(seq[p1:p1+k_size*order])


def mixedminstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int,
                    order: int = 2, strobe_fraction: float = 0.5) -> dict:
    """
    Mixed protocol to produce specified minstrobes and k-mers content fractions

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: a dictionary with positions along the string as keys and the mixed minstrobes/kmers as value
    """
    prime = 997
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order
    assert m_size + (order-1) * strobe_w_max_offset < len(seq), "Last minstrobes window position is exceeding the sequence length, consider to use a lower order"

    mixedminstrobes = {
        tuple(positions): h
        for positions, h in seq_to_mixedminstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w,
                order, denominator, numerator
            )
        }
    return mixedminstrobes


def mixedminstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, w: int, order: int = 2,
                         strobe_fraction: float = 1,
                         buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedminstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating minstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in mixedminstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order, strobe_fraction).items():
            yield p, m


def seq_to_hybridstrobes_iter(seq: str, k_size: int, w_min: int, w_max: int, w: int,
                              prime: int, order: int) -> Iterator[tuple]:
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
    hash_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)]
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
            for window_number, window in enumerate(tmp_index_dict[strobe_num]):
                # updating windows
                strobe_window, min_w, min_index, start, end = window
                new_w = hash_list[i + end]
                min_index, min_w = update_queue(
                    strobe_window, min_w, min_index, new_w, i, start, end
                )

                # update tmp_index_dict
                tmp_index_dict[strobe_num][window_number] = (
                    strobe_window, min_w, min_index, start, end
                )
                tmp_index.append((min_index, min_w))

            next_i, next_index_hash = tmp_index[index_hash % n_partition]
            positions.append(next_i)
            index_hash = index_hash + (strobe_num+1) * (-1)**(strobe_num+1) * next_index_hash

        yield positions, index_hash


def hybridstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                  strobe_w_max_offset: int, w: int, order: int = 2) -> dict:
    """
    Hybrid between minstrobes and randstrobes that uses both independent minimizers and a conditional dependence between strobes

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: a dictionary with positions along the string as keys and the hybridstrobes as value
    """
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"
    if k_size % order != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order

    if w == 1:
        hybridstrobes = {
            tuple(positions): index_hash
            for positions, index_hash in seq_to_hybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, prime, order
            )
        }

    else:  # thin out hybridstrobes
        hybridstrobes_tmp = [
            (positions, index_hash)
            for positions, index_hash in seq_to_hybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, prime, order
            )
        ]

        thinned_hybridstrobes = thinner(
            [index_hash for (positions, index_hash) in hybridstrobes_tmp],
            w
        )

        hybridstrobes = {}
        for p1, index_hash in thinned_hybridstrobes:
            if p1 < len(hybridstrobes_tmp):
                (positions, index_hash) = hybridstrobes_tmp[p1]
                hybridstrobes[tuple(positions)] = index_hash

    return hybridstrobes


def hybridstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                       strobe_w_max_offset: int, w: int, order: int = 2,
                       buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating hybridstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating hybridstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        # print(substring, len(substring))
        for p, m in hybridstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order).items():
            yield p, m


def seq_to_mixedhybridstrobes_iter(seq: str, k_size: int, w_min: int, w_max: int,
                                   w: int, prime: int, order: int, denominator: int,
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
    hash_list = [hash(seq[i:i+k_size]) for i in range(len(seq) - k_size + 1)]
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

    # print(tmp_index_dict)

    for i in range(len(hash_list) - w_max*order):  # temporary iteration
        index_hash = hash_list[i]
        positions = [i, ]

        for strobe_num in range(0, order-1):
            tmp_index = []
            for window_number, window in enumerate(tmp_index_dict[strobe_num]):
                # updating windows
                strobe_window, min_w, min_index, start, end = window
                new_w = hash_list[i + end]
                min_index, min_w = update_queue(
                    strobe_window, min_w, min_index, new_w, i, start, end
                )

                # update tmp_index_dict
                tmp_index_dict[strobe_num][window_number] = (
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


def mixedhybridstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                       strobe_w_max_offset: int, w: int,
                       order: int = 2, strobe_fraction: float = 0.5) -> dict:
    """
    Mixed protocol to produce specified hybridstrobes and k-mers content fractions

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: a dictionary with positions along the string as keys and the mixed hybridstrobes/kmers as value
    """
    prime = 997
    mixed_output = dict()
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order

    if w == 1:
        mixed_output = {
            tuple(positions): index_hash
            for positions, index_hash in seq_to_mixedhybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, prime,
                order, denominator, numerator
            )
        }

    else:  # thin out mixedhybridstrobes
        mixedhybridstrobes_tmp = [
            (positions, index_hash)
            for positions, index_hash in seq_to_mixedhybridstrobes_iter(
                seq, m_size, strobe_w_min_offset, strobe_w_max_offset, w, prime,
                order, denominator, numerator
                )
            ]

        thinned_mixedhybridstrobes = thinner(
            [index_hash for (positions, index_hash) in mixedhybridstrobes_tmp],
            w
        )

        hybridstrobes = {}
        for p1, index_hash in thinned_mixedhybridstrobes:
            if p1 < len(mixedhybridstrobes_tmp):
                (positions, index_hash) = mixedhybridstrobes_tmp[p1]
                mixed_output[tuple(positions)] = index_hash

    return mixed_output


def mixedhybridstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, w: int, order: int = 2,
                         strobe_fraction: float = 1,
                         buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedhybridstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating mixedhybridstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in mixedhybridstrobes(substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order, strobe_fraction).items():
            yield p, m


def seq_to_altstrobes_iter(seq: str, k_size1: int, k_size2: int, strobe_w_min_offset: int,
                           strobe_w_max_offset: int, prime: int, w: int,
                           order: int, generalized: bool=False, arg: int=50) -> Iterator[tuple]:
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
    assert order % 2 == 0, "Number of strobes has to be even in this implementation"

    hash_seq_list1 = [(i, hash(seq[i:i+k_size1])) for i in range(len(seq) - k_size1 + 1)]
    hash_seq_list2 = [(i, hash(seq[i:i+k_size2])) for i in range(len(seq) - k_size2 + 1)]
    strobe_w_min_offset1 = int(strobe_w_min_offset - (k_size1-k_size2)/2)
    strobe_w_max_offset1 = int(strobe_w_max_offset - (k_size1-k_size2)/2)
    strobe_w_min_offset2 = int(strobe_w_min_offset + (k_size1-k_size2)/2)
    strobe_w_max_offset2 = int(strobe_w_max_offset + (k_size1-k_size2)/2)
    required_length1 = order * strobe_w_max_offset - strobe_w_max_offset1
    required_length2 = order * strobe_w_max_offset - strobe_w_max_offset2

    # print(strobe_w_min_offset1, strobe_w_max_offset1)
    # print(strobe_w_min_offset2, strobe_w_max_offset2)

    # Number of distinct hash values
    # print(len(set(hash(seq[i:i+k_size1]) for i in range(len(seq) - k_size1 + 1))))
    # print(len(set(hash(seq[i:i+k_size2]) for i in range(len(seq) - k_size2 + 1))))

    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned1 = thinner([h for i, h in hash_seq_list1], w)
        hash_seq_list_thinned2 = thinner([h for i, h in hash_seq_list2], w)
    else:
        hash_seq_list_thinned1 = hash_seq_list1
        hash_seq_list_thinned2 = hash_seq_list2

    for (p1, hash_m1) in hash_seq_list_thinned1:  # [:-k_size]:
        if p1 >= int(len(hash_seq_list2) - order * (k_size1 + k_size2)/2 + k_size1):
            break
        # hash_m1 = hash_seq_list[p]
        # print(p1)

        if hash_m1 % 100 < arg: # decision about whether k1 or k2 should be seeded first
            seed_k1 = False
            index = [p1, ] # -p1
            min_hash_val = hash_m1
            required_length = required_length1
        else:
            seed_k1 = True
            if generalized:
                index = [p1, p1]
            else:
                index = [p1, p1+k_size1]
            hash_m1 = hash_seq_list2[p1][1]
            min_hash_val = hash_m1
            required_length = required_length2

        windows = list()
        if p1 + required_length <= min(len(hash_seq_list2), len(hash_seq_list1)):
            current_strobe_offset = p1
            for window_order in range(1, order):
                if seed_k1:  # seed strobe of length k1
                    start = strobe_w_min_offset1 + current_strobe_offset
                    current_strobe_offset += strobe_w_max_offset1  # update current offset
                    end = min(current_strobe_offset, len(hash_seq_list1))
                else:  # seed strobe of length k2
                    start = strobe_w_min_offset2 + current_strobe_offset
                    current_strobe_offset += strobe_w_max_offset2  # update current offset
                    end = min(current_strobe_offset, len(hash_seq_list2))

                windows.append((start, end))
                seed_k1 = not seed_k1  # alternate k1 and k2

        else:  # shrinking at end of sequence
            current_strobe_offset = p1
            current_strobe_length = p1
            for window_order in range(1, order):
                if seed_k1:
                    current_strobe_length += k_size1
                    start = (max(current_strobe_length,
                        len(hash_seq_list1) + strobe_w_min_offset1 - (order - window_order) * strobe_w_max_offset + int((k_size1-k_size2)/2)
                        )
                    )
                    current_strobe_offset += strobe_w_max_offset1  # update current offset
                    end = min(current_strobe_offset, len(hash_seq_list1))
                else:
                    current_strobe_length += k_size2
                    start = (max(
                        current_strobe_length,
                        len(hash_seq_list2) + strobe_w_min_offset2 - (order - window_order) * strobe_w_max_offset - int((k_size1-k_size2)/2)
                        )
                    )
                    current_strobe_offset += strobe_w_max_offset2  # update current offset
                    end = min(current_strobe_offset, len(hash_seq_list2))

                windows.append((start, end))
                seed_k1 = not seed_k1

        min_values = []
        for index_order in range(1, order):
            if len(range(*windows[index_order-1])) == 0:
                continue
            elif seed_k1:
                min_index, min_value = argmin([
                    (hash_m1 ^ hash_seq_list2[i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])
                min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list2[windows[index_order-1][0] + min_index][1]
                if generalized:
                    index.append(min_index+windows[index_order-1][0])
                    index.append(min_index+windows[index_order-1][0])
                else:
                    index.append(min_index+windows[index_order-1][0])
                    index.append(min_index+windows[index_order-1][0]+k_size1)
                min_values.append(min_value)
            else:
                min_index, min_value = argmin([
                    (hash_m1 ^ hash_seq_list1[i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])
                min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list1[windows[index_order-1][0] + min_index][1]
                index.append(min_index+windows[index_order-1][0]) # -
                min_values.append(min_value)
            seed_k1 = not seed_k1

        yield index, min_hash_val #, min_values


def altstrobes(seq: str, k_size: int, strobe_w_min_offset: int,
               strobe_w_max_offset: int, w: int, order: int = 2,
               prime: int = 997, arg: int = 50) -> dict:
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
    assert order % 2 == 0, "Number of strobes has to be even in this implementation"

    if k_size % int(1.5*order) != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, int(order*1.5), k_size - k_size % int(1.5*order)))
        k_size = k_size - k_size % int(1.5*order)
    m_size1 = int(2*k_size/(3*order))
    m_size2 = int(4*k_size/(3*order))

    altstrobes = {tuple(index): h for index, h in seq_to_altstrobes_iter(
        seq, m_size1, m_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, arg=arg
    )}
    return altstrobes


def altstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int, order: int = 2,
                    buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating altstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating altstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in altstrobes(
                substring, k_size, strobe_w_min_offset, strobe_w_max_offset,
                w, order).items():

            yield p, m


def altstrobes_generalized(seq: str, k_size1: int, k_size2: int, strobe_w_min_offset: int,
               strobe_w_max_offset: int, w: int, order: int = 2,
               prime: int = 997, arg: int = 50) -> dict:
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
    assert order == 2, "Altstrobes are just implemented for order 2"

    altstrobes = {tuple(index): h for index, h in seq_to_altstrobes_iter(
        seq, k_size1, k_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, generalized=True, arg=arg
    )}
    return altstrobes


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
    assert order % 2 == 0, "Number of strobes has to be even in this implementation"

    hash_seq_list1 = [(i, hash(seq[i:i+k_size1])) for i in range(len(seq) - k_size1 + 1)]
    hash_seq_list2 = [(i, hash(seq[i:i+k_size2])) for i in range(len(seq) - k_size2 + 1)]
    strobe_w_min_offset1 = int(strobe_w_min_offset - (k_size1-k_size2)/2)
    strobe_w_max_offset1 = int(strobe_w_max_offset - (k_size1-k_size2)/2)
    strobe_w_min_offset2 = int(strobe_w_min_offset + (k_size1-k_size2)/2)
    strobe_w_max_offset2 = int(strobe_w_max_offset + (k_size1-k_size2)/2)
    required_length1 = order * strobe_w_max_offset - strobe_w_max_offset1
    required_length2 = order * strobe_w_max_offset - strobe_w_max_offset2

    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        hash_seq_list_thinned1 = thinner([h for i, h in hash_seq_list1], w)
        hash_seq_list_thinned2 = thinner([h for i, h in hash_seq_list2], w)
    else:
        hash_seq_list_thinned1 = hash_seq_list1
        hash_seq_list_thinned2 = hash_seq_list2

    for (p1, hash_m1) in hash_seq_list_thinned1:  # [:-k_size]:
        if p1 >= int(len(hash_seq_list2) - order * (k_size1 + k_size2)/2 + k_size1):
            break
        # hash_m1 = hash_seq_list[p]
        # print(p1)

        if (hash_m1//100) % denominator < numerator:  # pick altstrobe
            if hash_m1 % 2 == 0: # decision about whether k1 or k2 should be seeded first
                seed_k1 = False
                index = [p1, ] # -p1
                min_hash_val = hash_m1
                required_length = required_length1
            else:
                seed_k1 = True
                index = [p1, p1+k_size1]
                hash_m1 = hash_seq_list2[p1][1]
                min_hash_val = hash_m1
                required_length = required_length2

            windows = list()
            if p1 + required_length <= min(len(hash_seq_list2), len(hash_seq_list1)):
                current_strobe_offset = p1
                for window_order in range(1, order):
                    if seed_k1:  # seed strobe of length k1
                        start = strobe_w_min_offset1 + current_strobe_offset
                        current_strobe_offset += strobe_w_max_offset1  # update current offset
                        end = min(current_strobe_offset, len(hash_seq_list1))
                    else:  # seed strobe of length k2
                        start = strobe_w_min_offset2 + current_strobe_offset
                        current_strobe_offset += strobe_w_max_offset2  # update current offset
                        end = min(current_strobe_offset, len(hash_seq_list2))

                    windows.append((start, end))
                    seed_k1 = not seed_k1  # alternate k1 and k2

            else:  # shrinking at end of sequence
                current_strobe_offset = p1
                current_strobe_length = p1
                for window_order in range(1, order):
                    if seed_k1:
                        current_strobe_length += k_size1
                        start = (max(current_strobe_length,
                            len(hash_seq_list1) + strobe_w_min_offset1 - (order - window_order) * strobe_w_max_offset + int((k_size1-k_size2)/2)
                            )
                        )
                        current_strobe_offset += strobe_w_max_offset1  # update current offset
                        end = min(current_strobe_offset, len(hash_seq_list1))
                    else:
                        current_strobe_length += k_size2
                        start = (max(
                            current_strobe_length,
                            len(hash_seq_list2) + strobe_w_min_offset2 - (order - window_order) * strobe_w_max_offset - int((k_size1-k_size2)/2)
                            )
                        )
                        current_strobe_offset += strobe_w_max_offset2  # update current offset
                        end = min(current_strobe_offset, len(hash_seq_list2))

                    windows.append((start, end))
                    seed_k1 = not seed_k1

            min_values = []
            for index_order in range(1, order):
                if seed_k1:
                    min_index, min_value = argmin([
                        (hash_m1 ^ hash_seq_list2[i][1]) % prime
                        for i in range(*windows[index_order-1])
                    ])
                    min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list2[windows[index_order-1][0] + min_index][1]
                    index.append(min_index+windows[index_order-1][0])
                    index.append(min_index+windows[index_order-1][0]+k_size1)
                    min_values.append(min_value)
                else:
                    min_index, min_value = argmin([
                        (hash_m1 ^ hash_seq_list1[i][1]) % prime
                        for i in range(*windows[index_order-1])
                    ])
                    min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list1[windows[index_order-1][0] + min_index][1]
                    index.append(min_index+windows[index_order-1][0]) # -
                    min_values.append(min_value)
                seed_k1 = not seed_k1

            yield index, min_hash_val #, min_values

        else:  # pick k-mer
            index = tuple(p1 + (strobe_num) * k_size1 for strobe_num in range(int(1.5*order)))
            l = int((k_size1+k_size2)/2*order)
            yield index, hash(seq[p1: p1+l])


def altstrobes_generalized_iter(seq: str, k_size1: int, k_size2: int, strobe_w_min_offset: int,
                    strobe_w_max_offset: int, w: int, order: int = 2,
                    buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating altstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating altstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in altstrobes_generalized(
                substring, k_size1, k_size2, strobe_w_min_offset, strobe_w_max_offset,
                w, order=order).items():

            yield p, m


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
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, int(order*1.5), k_size - k_size % int(1.5*order)))
        k_size = k_size - k_size % int(1.5*order)

    m_size1 = int(2*k_size/(3*order))
    m_size2 = int(4*k_size/(3*order))

    mixedaltstrobes = {tuple(index): h for index, h in seq_to_mixedaltstrobes_iter(
        seq, m_size1, m_size2, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator
    )}

    return mixedaltstrobes


def mixedaltstrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                     strobe_w_max_offset: int, w: int, order: int = 2, strobe_fraction: float = 1,
                     buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedaltstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating mixedaltstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in mixedaltstrobes(
                substring, k_size, strobe_w_min_offset, strobe_w_max_offset,
                w, order, strobe_fraction).items():

            yield p, m


def seq_to_multistrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                                     strobe_w_max_offset: int, prime: int, w: int,
                                     order: int, k_boundary: int,
                                     arg: int=50) -> Iterator[tuple]:
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
    assert k_boundary > 0, "[Error] It is not possible to sample 0mers"
    assert order % 2 == 0, "Number of strobes has to be even in this implementation"

    k_size = int(2*k_size/order)

    hash_seq_lists = []
    strobe_w_min_offsets = []
    strobe_w_max_offsets = []
    hash_seq_lists_thinned = []
    required_length = (order-1) * strobe_w_max_offset + k_size/2

    for k1 in range(k_boundary, k_size-k_boundary+1):
        hash_seq_lists.append([(i, hash(seq[i:i+k1])) for i in range(len(seq) - k1 + 1)])
        strobe_w_min_offsets.append(int(strobe_w_min_offset - k1 + k_size/2))
        strobe_w_max_offsets.append(int(strobe_w_max_offset - k1 + k_size/2))

    # thinning
    if w > 1:
        # produce a subset of positions, still with same index as in full sequence
        for n in range(len(hash_seq_lists)):
            hash_seq_lists_thinned.append(thinner([h for i, h in hash_seq_lists[n]], w))
    else:
        for n in range(len(hash_seq_lists)):
            hash_seq_lists_thinned.append(hash_seq_lists[n])

    if len(hash_seq_lists) % 2 == 0:
        k_size_options = int(len(hash_seq_lists)/2)
    else:
        k_size_options = int(len(hash_seq_lists)/2+1)

    for p1, (p1, hash_m1) in enumerate(hash_seq_lists_thinned[0]):  # [:-k_size]:
        if p1 >= min(len(hash_seq_lists[0]) - (order-1)*k_boundary, len(hash_seq_lists[-1]) - (order-1)*(k_size-k_boundary)):
            break
        # hash_m1 = hash_seq_list[p]

        k_size_selection = hash_m1 % k_size_options

        if hash_seq_lists[k_size_selection][p1][1] // 100 % 100 < arg:
            k_size_2_selection = -(k_size_selection+1)
            k_size1 = k_size_selection + k_boundary
            k_size2 = k_size - k_size1
        else:
            k_size_2_selection = hash_m1 % k_size_options
            k_size_selection = -(k_size_2_selection+1)
            k_size2 = k_size_2_selection + k_boundary
            k_size1 = k_size - k_size2
        hash_m1 = hash_seq_lists[k_size_selection][p1][1]

        if p1 + required_length <= len(seq):
            windows = list()
            current_strobe_offset = p1
            for window_order in range(1, order):
                if window_order % 2 == 1:
                    start = strobe_w_min_offsets[k_size_2_selection] + current_strobe_offset
                    current_strobe_offset += strobe_w_max_offsets[k_size_2_selection]
                    end = min(current_strobe_offset, len(hash_seq_lists[k_size_2_selection]))
                else:
                    start = strobe_w_min_offsets[k_size_selection] + current_strobe_offset
                    current_strobe_offset += strobe_w_max_offsets[k_size_selection]
                    end = min(current_strobe_offset, len(hash_seq_lists[k_size_selection]))
                windows.append((start, end))

        else:
            windows = list()
            current_strobe_offset = p1
            current_strobe_length = p1
            for window_order in range(1, order):
                if window_order % 2 == 1:
                    current_strobe_length += k_size1
                    start = (max(
                        current_strobe_length,
                        len(hash_seq_lists[k_size_2_selection]) + strobe_w_min_offsets[k_size_2_selection] - (order - window_order) * strobe_w_max_offset
                        )
                    )
                    current_strobe_offset += strobe_w_max_offsets[k_size_2_selection]
                    end = min(current_strobe_offset, len(hash_seq_lists[k_size_2_selection]))
                else:
                    current_strobe_length += k_size2
                    start = (max(
                        current_strobe_length,
                        len(hash_seq_lists[k_size_selection]) + strobe_w_min_offsets[k_size_selection] - (order - window_order) * strobe_w_max_offset
                        )
                    )
                    current_strobe_offset += strobe_w_max_offsets[k_size_selection]
                    end = min(current_strobe_offset, len(hash_seq_lists[k_size_selection]))

                windows.append((start, end))

        index = [range(p1, p1+k_size1), ] # -p1
        min_values = []
        min_hash_val = hash_m1
        for index_order in range(1, order):
            if len(range(*windows[index_order-1])) == 0:
                continue
            elif index_order % 2 == 1:
                min_index, min_value = argmin([
                    (hash_m1 ^ hash_seq_lists[k_size_2_selection][i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])

                min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_lists[k_size_2_selection][windows[index_order-1][0] + min_index][1]
                index.append(range(min_index+windows[index_order-1][0], min_index+windows[index_order-1][0]+k_size2))
                min_values.append(min_value)
            else:
                min_index, min_value = argmin([
                    (hash_m1 ^ hash_seq_lists[k_size_selection][i][1]) % prime
                    for i in range(*windows[index_order-1])
                ])

                min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_lists[k_size_selection][windows[index_order-1][0] + min_index][1]
                index.append(range(min_index+windows[index_order-1][0], min_index+windows[index_order-1][0]+k_size1))
                min_values.append(min_value)

        if min_hash_val == 0:
            print(len(seq), index, start, end, k_size1, k_size2, min_hash_val)
        yield index, min_hash_val# , min_values


def multistrobes(seq: str, k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, w: int, order: int = 2,
                         prime: int = 997, k_boundary: int = 5,
                         arg: int = 50) -> dict:
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

    multistrobes = {tuple(index): h for index, h in seq_to_multistrobes_iter(
        seq, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, k_boundary, arg
    )}
    return multistrobes


def multistrobes_iter(seq: str, k_size: int, strobe_w_min_offset: int,
                              strobe_w_max_offset: int, w: int, order: int = 2,
                              k_boundary: int=5, buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating altstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating altstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i:i+buffer_size]
        for p, m in multistrobes(
                substring, k_size, strobe_w_min_offset, strobe_w_max_offset,
                w, order=order, k_boundary=k_boundary).items():

            yield p, m


"""
Mixedstrobemer techniques sampling seeds using two strobemer protocols
"""

def seq_to_mixedstrobes_iter(method1: str, method2: str, seq: str,
                             k_size: int, strobe_w_min_offset: int,
                             strobe_w_max_offset: int, prime: int, w: int,
                             order: int, denominator: int,
                             numerator: int) -> Iterator[tuple]:
    """
    Iterator for creating of mixedmixed of any order

    :param seq: a string with a nucleotide sequence
    :param k_size: length of each strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param denominator: denominator and numerator determine the fraction of sampled strobemers
    :param numerator: denominator and numerator determine the fraction of sampled strobemers
    :returns: an iterator for creating mixedstrobes
    """
    hash_seq_list = [(i, hash(seq[i:i+k_size])) for i in range(len(seq) - k_size + 1)]
    # thinning
    if w > 1:
        # produce a subset of positions, still with samme index as in full sequence
        hash_seq_list_thinned = thinner([h for i, h in hash_seq_list], w)
    else:
        hash_seq_list_thinned = hash_seq_list

    if "minstrobes" in (method1, method2):
        strobes = deque(thinner([h for i, h in hash_seq_list], strobe_w_max_offset - strobe_w_min_offset))
        strobes_dict = {strobe_num: copy.deepcopy(strobes) for strobe_num in range(1, order)}
    if "hybridstrobes" in (method1, method2):
        n_partition = 3
        w_p = (strobe_w_max_offset - strobe_w_min_offset) // n_partition
        tmp_index_dict = dict()
        for strobe_num in range(0, order-1):
            tmp_index = []
            for partition in range(0, n_partition):
                start = strobe_w_min_offset + strobe_w_min_offset*strobe_num + w_p*partition
                end = (
                    strobe_w_max_offset + strobe_w_max_offset*strobe_num if partition + 1 == n_partition
                    else strobe_w_min_offset + strobe_w_max_offset*strobe_num + w_p + w_p*partition
                )
                strobe_window = deque(list(chain(zip(*hash_seq_list[start: end])))[1])
                min_index, min_w = argmin(strobe_window)
                min_index = min_index + strobe_w_min_offset + strobe_w_max_offset*strobe_num + w_p*partition
                tmp_index.append(
                    (
                        strobe_window,
                        min_w,
                        min_index,
                        start, end
                    )
                )
            tmp_index_dict[strobe_num] = tmp_index
        # print(tmp_index_dict)
    if "altstrobes" in (method1, method2):
        k_size1 = int(k_size * 2 / 3)
        k_size2 = int(k_size * 4 / 3)
        hash_seq_list1 = [(i, hash(seq[i:i+k_size1])) for i in range(len(seq) - k_size1 + 1)]
        hash_seq_list2 = [(i, hash(seq[i:i+k_size2])) for i in range(len(seq) - k_size2 + 1)]
        strobe_w_min_offset1 = int(strobe_w_min_offset - (k_size1-k_size2)/2)
        strobe_w_max_offset1 = int(strobe_w_max_offset - (k_size1-k_size2)/2)
        strobe_w_min_offset2 = int(strobe_w_min_offset + (k_size1-k_size2)/2)
        strobe_w_max_offset2 = int(strobe_w_max_offset + (k_size1-k_size2)/2)
        required_length1 = order * strobe_w_max_offset - strobe_w_max_offset1
        required_length2 = order * strobe_w_max_offset - strobe_w_max_offset2
    if "multistrobes" in (method1, method2):
        k_boundary = 5
        hash_seq_lists = []
        strobe_w_min_offsets = []
        strobe_w_max_offsets = []
        for k1 in range(k_boundary, order*k_size-k_boundary+1):
            hash_seq_lists.append([(i, hash(seq[i:i+k1])) for i in range(len(seq) - k1 + 1)])
            strobe_w_min_offsets.append(int(strobe_w_min_offset + k1 - 15))
            strobe_w_max_offsets.append(int(strobe_w_max_offset + k1 - 15))
        if len(hash_seq_lists) % 2 == 0:
            k_size_options = int(len(hash_seq_lists)/2)
        else:
            k_size_options = int(len(hash_seq_lists)/2+1)

    for (p1, hash_m1) in hash_seq_list_thinned:  # [:-k_size]:
        if p1 >= len(hash_seq_list) - (order-1)*k_size:
            break

        if hash_m1//100 % denominator < numerator:  # pick method1
            if method1 == "randstrobes":
                index, min_hash_val = sample_randstrobes(seq, hash_seq_list, p1, hash_m1, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, order)
            elif method1 == "minstrobes":
                index, min_hash_val = sample_minstrobes(seq, hash_seq_list, p1, hash_m1, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, order, strobes_dict)
            elif method1 == "hybridstrobes":
                if p1 >= len(hash_seq_list) - strobe_w_max_offset*order:
                    break
                index, min_hash_val, tmp_index_dict = sample_hybridstrobes(seq, hash_seq_list, p1, hash_m1, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, order, n_partition, tmp_index_dict)
            elif method1 == "altstrobes":
                if p1 >= len(hash_seq_list2) - (order-1)*k_size2:
                    break
                index, min_hash_val = sample_altstrobes(seq, hash_seq_list1, hash_seq_list2, p1, hash_m1, k_size1, k_size2, strobe_w_min_offset1, strobe_w_max_offset1,
                                                        strobe_w_min_offset2, strobe_w_max_offset2, strobe_w_max_offset, prime, order, required_length1, required_length2)
            elif method1 == "multistrobes":
                if p1 >= min(len(hash_seq_lists[0]) - (order-1)*k_boundary, len(hash_seq_lists[-1]) - (order-1)*(order*k_size-k_boundary)):
                    break
                index, min_hash_val = sample_multistrobes(seq, hash_seq_lists, p1, hash_m1, order*k_size, k_size_options, strobe_w_min_offsets, strobe_w_max_offsets, prime, order, k_boundary)
            else:
                print("Seeding Technique not known")
                raise NotImplementedError

            yield index, min_hash_val

        else:  # pick method2
            if method2 == "randstrobes":
                index, min_hash_val = sample_randstrobes(seq, hash_seq_list, p1, hash_m1, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, order)
            elif method2 == "minstrobes":
                index, min_hash_val = sample_minstrobes(seq, hash_seq_list, p1, hash_m1, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, order, strobes_dict)
            elif method2 == "hybridstrobes":
                if p1 >= len(hash_seq_list) - strobe_w_max_offset*order:
                    break
                index, min_hash_val, tmp_index_dict = sample_hybridstrobes(seq, hash_seq_list, p1, hash_m1, k_size, strobe_w_min_offset, strobe_w_max_offset, prime, order, n_partition, tmp_index_dict)
            elif method2 == "altstrobes":
                if p1 >= len(hash_seq_list2) - (order-1)*k_size2:
                    break
                index, min_hash_val = sample_altstrobes(seq, hash_seq_list1, hash_seq_list2, p1, hash_m1, k_size1, k_size2, strobe_w_min_offset1, strobe_w_max_offset1,
                                                        strobe_w_min_offset2, strobe_w_max_offset2, strobe_w_max_offset, prime, order, required_length1, required_length2)
            elif method2 == "multistrobes":
                if p1 >= min(len(hash_seq_lists[0]) - (order-1)*k_boundary, len(hash_seq_lists[-1]) - (order-1)*(order*k_size-k_boundary)):
                    break
                index, min_hash_val = index, min_hash_val = sample_multistrobes(seq, hash_seq_lists, p1, hash_m1, order*k_size, k_size_options, strobe_w_min_offsets, strobe_w_max_offsets, prime, order, k_boundary)
            else:
                print("Seeding Technique not known")
                raise NotImplementedError

            yield index, min_hash_val


def sample_randstrobes(seq: str, hash_seq_list: list, p1: int, hash_m1: int,
                       k_size: int, strobe_w_min_offset: int, strobe_w_max_offset: int,
                       prime: int, order: int):
    """
    Sample one randstrobe at a defined position

    :param hash_seq_list: list with tuples of positions and corresponding hash values
    :param p1: start position of first strobe
    :param hash_m1: hash of first strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: one sampled randstrobe
    """
    if p1 + (order-1) * strobe_w_max_offset <= len(hash_seq_list):
        windows = list()
        for window_order in range(1, order):
            start = p1 + strobe_w_min_offset + (window_order-1) * strobe_w_max_offset
            end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
            windows.append((start, end))

    else:
        windows = list()
        for window_order in range(1, order):
            start = (max(
                p1+window_order*k_size,
                len(hash_seq_list) + strobe_w_min_offset - (order - window_order) * strobe_w_max_offset
                )
            )

            end = min(p1 + window_order * strobe_w_max_offset, len(hash_seq_list))
            windows.append((start, end))

    index = [range(p1, p1+k_size), ]
    min_hash_val = hash_m1
    for index_order in range(1, order):
        min_index, min_value = argmin([
            (hash_m1 ^ hash_seq_list[i][1])
            for i in range(*windows[index_order-1])
        ])

        min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list[windows[index_order-1][0] + min_index][1]
        index.append(range(min_index+windows[index_order-1][0], min_index+windows[index_order-1][0]+k_size))
    return index, min_hash_val


def sample_hybridstrobes(seq: str, hash_seq_list: list, p1: int, hash_m1: int,
                         k_size: int, strobe_w_min_offset: int,
                         strobe_w_max_offset: int, prime: int, order: int,
                         n_partition: int, tmp_index_dict: dict):
    """
    Sample one hybridstrobe at a defined position

    :param hash_seq_list: list with tuples of positions and corresponding hash values
    :param p1: start position of first strobe
    :param hash_m1: hash of first strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param n_partition: number of partitions
    :returns: one sampled hybridstrobe
    """
    positions = [range(p1, p1+k_size), ]
    hash_value = hash_m1

    for strobe_num in range(0, order-1):
        tmp_index = []
        for window_number, window in enumerate(tmp_index_dict[strobe_num]):
            # updating windows
            strobe_window, min_w, min_index, start, end = window
            new_w = hash_seq_list[p1 + end][1]
            min_index, min_w = update_queue(
                strobe_window, min_w, min_index, new_w, p1, start, end
            )

            # update tmp_index_dict
            tmp_index_dict[strobe_num][window_number] = (
                strobe_window, min_w, min_index, start, end
            )
            tmp_index.append((min_index, min_w))

        next_i, next_index_hash = tmp_index[hash_value % n_partition]
        positions.append(range(next_i, next_i+k_size))
        hash_value += (strobe_num+1) * (-1)**(strobe_num+1) * next_index_hash

    return positions, hash_value, tmp_index_dict


def sample_minstrobes(seq: str, hash_seq_list: list, p1: int, hash_m1: int,
                      k_size: int, strobe_w_min_offset: int, strobe_w_max_offset: int,
                      prime: int, order: int, strobes_dict: dict):
    """
    Sample one minstrobe at a defined position

    :param hash_seq_list: list with tuples of positions and corresponding hash values
    :param p1: start position of first strobe
    :param hash_m1: hash of first strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobes_dict:
    :returns: one sampled minstrobe
    """
    positions = [range(p1, p1+k_size), ]
    hash_value = hash_m1

    for strobe_num in range(1, order):
        if p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset < len(seq):
            while strobes_dict[strobe_num][0][0] < min(p1 + k_size + strobe_w_min_offset + (strobe_num-1) * strobe_w_max_offset, len(hash_seq_list)-1):
                l = strobes_dict[strobe_num].popleft()
        p, hash_val = strobes_dict[strobe_num][0]
        positions.append(range(p, p+k_size))
        hash_value += strobe_num * (-1)**strobe_num * hash_val
    return positions, hash_value


def sample_altstrobes(seq: str, hash_seq_list1: list, hash_seq_list2: list, p1: int,
                      hash_m1: int, k_size1: int, k_size2: int, strobe_w_min_offset1: int,
                      strobe_w_max_offset1: int, strobe_w_min_offset2: int,
                      strobe_w_max_offset2: int, strobe_w_max_offset: int, prime: int, order: int,
                      required_length1: int, required_length2: int):
    """
    Sample one altstrobe at a defined position

    :param hash_seq_list: list with tuples of positions and corresponding hash values
    :param p1: start position of first strobe
    :param hash_m1: hash of first strobe
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param prime: prime number (q) in minimizing h(m)+h(mj) mod q
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :returns: one sampled altstrobe
    """
    hash_m1 = hash_seq_list1[p1][1]

    if hash_m1 % 2 == 0: # first x, then 2x (e.g. 10-20)
        seed_k1 = False
        index = [range(p1, p1+k_size1), ] # -p1
        min_hash_val = hash_m1
        required_length = required_length1
    else:
        seed_k1 = True
        index = [range(p1, p1+k_size2), ]
        hash_m1 = hash_seq_list2[p1][1]
        min_hash_val = hash_m1
        required_length = required_length2

    windows = list()
    if p1 + required_length <= min(len(hash_seq_list2), len(hash_seq_list1)):
        current_strobe_offset = p1
        for window_order in range(1, order):
            if seed_k1:  # seed strobe of length k1
                start = strobe_w_min_offset1 + current_strobe_offset
                current_strobe_offset += strobe_w_max_offset1  # update current offset
                end = min(current_strobe_offset, len(hash_seq_list1))
            else:  # seed strobe of length k2
                start = strobe_w_min_offset2 + current_strobe_offset
                current_strobe_offset += strobe_w_max_offset2  # update current offset
                end = min(current_strobe_offset, len(hash_seq_list2))

            windows.append((start, end))
            seed_k1 = not seed_k1  # alternate k1 and k2

    else:  # shrinking at end of sequence
        current_strobe_offset = p1
        current_strobe_length = p1
        for window_order in range(1, order):
            if seed_k1:
                current_strobe_length += k_size1
                start = (max(current_strobe_length,
                    len(hash_seq_list1) + strobe_w_min_offset1 - (order - window_order) * strobe_w_max_offset + int((k_size1-k_size2)/2)
                    )
                )
                current_strobe_offset += strobe_w_max_offset1  # update current offset
                end = min(current_strobe_offset, len(hash_seq_list1))
            else:
                current_strobe_length += k_size2
                start = (max(
                    current_strobe_length,
                    len(hash_seq_list2) + strobe_w_min_offset2 - (order - window_order) * strobe_w_max_offset - int((k_size1-k_size2)/2)
                    )
                )
                current_strobe_offset += strobe_w_max_offset2  # update current offset
                end = min(current_strobe_offset, len(hash_seq_list2))

            windows.append((start, end))
            seed_k1 = not seed_k1

    min_values = []
    for index_order in range(1, order):
        if seed_k1:
            min_index, min_value = argmin([
                (hash_m1 ^ hash_seq_list2[i][1]) % prime
                for i in range(*windows[index_order-1])
            ])
            min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list2[windows[index_order-1][0] + min_index][1]
            index.append(range(min_index+windows[index_order-1][0], min_index+windows[index_order-1][0]+k_size2))
            # min_values.append(min_value)
        else:
            min_index, min_value = argmin([
                (hash_m1 ^ hash_seq_list1[i][1]) % prime
                for i in range(*windows[index_order-1])
            ])
            min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_list1[windows[index_order-1][0] + min_index][1]
            index.append(range(min_index+windows[index_order-1][0], min_index+windows[index_order-1][0]+k_size1)) # -
            # min_values.append(min_value)

    return index, min_hash_val


def sample_multistrobes(seq: str, hash_seq_lists: list, p1: int,
                        hash_m1: int, k_size: int, k_size_options: int,
                        strobe_w_min_offsets: list,
                        strobe_w_max_offsets: list, prime: int,
                        order: int, k_boundary: int):
    """
    """
    hash_m1 = hash_seq_lists[0][p1][1]  # not necessary, but guarantees identical values to multistrobes
    k_size_selection = hash_m1 % k_size_options
    k_size_2_selection = -(k_size_selection+1)
    k_size1 = k_size_selection + k_boundary
    k_size2 = k_size - k_size1
    hash_m1 = hash_seq_lists[k_size_selection][p1][1]
    # print(hash_m1, k_size_options, k_size1, k_size2)

    if hash_m1 // 100 % 2 == 0: # first k_size1, then k_size2 (e.g. 10-20)
        # print("10-20")
        if p1 + (order-1) * strobe_w_max_offsets[k_size_2_selection] <= len(hash_seq_lists[k_size_2_selection]):
            windows = list()
            for window_order in range(1, order):
                start = p1 + strobe_w_min_offsets[k_size_2_selection] + (window_order-1) * strobe_w_max_offsets[k_size_2_selection]
                end = min(p1 + window_order * strobe_w_max_offsets[k_size_2_selection], len(hash_seq_lists[k_size_2_selection]))
                windows.append((start, end))

        else:
            windows = list()
            for window_order in range(1, order):
                start = (max(
                    p1+window_order*k_size2,
                    len(hash_seq_lists[k_size_2_selection]) + strobe_w_min_offsets[k_size_2_selection] - (order - window_order) * strobe_w_max_offsets[k_size_2_selection]
                    )
                )

                end = min(p1 + window_order * strobe_w_max_offsets[k_size_2_selection], len(hash_seq_lists[k_size_2_selection]))
                windows.append((start, end))

        index = [range(p1, p1+k_size1), ] # -p1
        min_values = []
        min_hash_val = hash_m1
        for index_order in range(1, order):
            min_index, min_value = argmin([
                (hash_m1 ^ hash_seq_lists[k_size_2_selection][i][1]) % prime
                for i in range(*windows[index_order-1])
            ])

            min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_lists[k_size_2_selection][windows[index_order-1][0] + min_index][1]
            index.append(range(min_index+windows[index_order-1][0], min_index+windows[index_order-1][0]+k_size2))
            min_values.append(min_value)

    else: # first k_size2, then k_size1 (e.g. 20-10)
        # print("20-10")

        if p1 + (order-1) * strobe_w_max_offsets[k_size_selection] <= len(hash_seq_lists[k_size_selection]):
            windows = list()
            for window_order in range(1, order):
                start = p1 + strobe_w_min_offsets[k_size_selection] + (window_order-1) * strobe_w_max_offsets[k_size_selection]
                end = min(p1 + window_order * strobe_w_max_offsets[k_size_selection], len(hash_seq_lists[k_size_selection]))
                windows.append((start, end))

        else:
            windows = list()
            for window_order in range(1, order):
                start = (max(
                    p1+window_order*k_size1,
                    len(hash_seq_lists[k_size_selection]) + strobe_w_min_offsets[k_size_selection] - (order - window_order) * strobe_w_max_offsets[k_size_selection]
                    )
                )

                end = min(p1 + window_order * strobe_w_max_offsets[k_size_selection], len(hash_seq_lists[k_size_selection]))
                windows.append((start, end))

        index = [range(p1, p1+k_size2)]
        min_values = []
        hash_m1 = hash_seq_lists[k_size_2_selection][p1][1]
        min_hash_val = hash_m1
        for index_order in range(1, order):
            min_index, min_value = argmin([
                (hash_m1 ^ hash_seq_lists[k_size_selection][i][1]) % prime
                for i in range(*windows[index_order-1])
            ])

            min_hash_val = min_hash_val + ((index_order+1) * (-1)**index_order) * hash_seq_lists[k_size_selection][windows[index_order-1][0] + min_index][1]
            index.append(range(min_index+windows[index_order-1][0], min_index+windows[index_order-1][0]+k_size1)) # -
            min_values.append(min_value)

    # print(index, start, end, min_hash_val)
    return index, min_hash_val


def mixedstrobes(method1: str, method2: str, seq: str, k_size: int,
                 strobe_w_min_offset: int, strobe_w_max_offset: int, w: int,
                 order: int = 2, strobe_fraction: float = 0.5) -> dict:
    """
    Mixed protocol to produce specified randstrobes and k-mers content fractions

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param strobe_fraction: a fraction of sampled strobemers, rest kmers (0 <= strobe_fraction <= 1)
    :returns: a dictionary with positions along the string as keys and the mixed randstrobes/kmers as value
    """
    fraction = Fraction(str(strobe_fraction))
    denominator = fraction.denominator
    numerator = fraction.numerator
    prime = 997
    assert strobe_w_min_offset > 0, "Minimum strobemer offset has to be greater than 0 in this implementation"

    if k_size % order != 0:
        print("WARNING: kmer size ({0}) is not evenly divisible with {1}, will use {2} as kmer size: ".format(k_size, order, k_size - k_size % order))
        k_size = k_size - k_size % order
    m_size = k_size//order

    mixedstrobes = {
        tuple(index): h
        for index, h in seq_to_mixedstrobes_iter(
            method1, method2, seq, m_size, strobe_w_min_offset, strobe_w_max_offset, prime, w, order, denominator, numerator
        )
    }
    return mixedstrobes


def mixedstrobes_iter(method1: str, method2: str, seq: str, k_size: int,
                      strobe_w_min_offset: int, strobe_w_max_offset: int,
                      w: int, order: int = 2, strobe_fraction: float = 1,
                      buffer_size: int = 10000000) -> Iterator[tuple]:
    """
    Generator for creating mixedstrobes (less memory requiring)

    :param seq: a string with a nucleotide sequence
    :param k_size: length of all strobes (len(strobe_1) +  ... + len(strobe_n))
    :param strobe_w_min_offset: minimum window offset to the previous window (wMin > 0)
    :param strobe_w_max_offset: maximum window offset to the previous window (wMin <= wMax)
    :param w: number of hashes used in a sliding window for thinning (w=1 means no thinning)
    :param order: number of substrings/strobes
    :param buffer_size: size of buffer at the end of the sequence
    :returns: an iterator for creating mixedstrobes
    """
    for i in range(0, len(seq), buffer_size):
        substring = seq[i: i+buffer_size]
        for p, m in mixedstrobes(method1, method2, substring, k_size, strobe_w_min_offset, strobe_w_max_offset, w, order, strobe_fraction).items():
            yield p, m
