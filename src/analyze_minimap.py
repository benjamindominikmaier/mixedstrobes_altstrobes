#!/usr/bin/env python3.9
# -*- coding: utf-8 -*

"""Compute fraction of correctly mapped reads from minimap2 results

This script takes a SAM output file from minimap2-strobemers (--input_file) and
computes the fraction of correctly mapped reads (as defined by the existence of
an overlap between query and reference position and right direction) and the
average (mean) number of correctly mapped nucleotides.
"""

__authors__ = ["Benjamin D. Maier"]
__copyright__ = "Copyright Benjamin D. Maier & Kristoffer Sahlin | Sahlin Group"
__organization__ = "Department of Mathematics, Science for Life Laboratory, Stockholm University, 106 91, Stockholm, Sweden."
__credits__ = ["Benjamin D. Maier & Kristoffer Sahlin"]
__contact__ = "bmaier [at] ebi.ac.uk"
__date__ = "2023/03/10"
__created__ = "2022/07/25"
__deprecated__ = False
__license__ = "MIT"
__maintainer__ = "Kristoffer Sahlin"
__email__ = "kristoffer.sahlin [at] scilifelab.se"
__status__ = "DSML Lvl. 1 - Concept"

import os
import sys
import csv
import argparse
import random
from itertools import chain
import more_itertools as mit


def main(args):
    output_dict = {}
    counter_fw = 0
    counter_rev = 0
    with open(args.input_file, newline = "\n") as SAM:
        for line in csv.reader(SAM, delimiter="\t"):
            if not line[0].startswith("@"):
                # print(line[0:3])
                if line[1] == "0":  # primary alignment fw
                    align_length = len(line[9])
                    nr, chr, dir, q_pos, mut = line[0].split("|")
                    # print(nr, chr, dir, q_pos, mut)
                    q_pos = int(q_pos)
                    r_pos = int(line[3])
                    if (dir=="forw") & (line[2] == chr) & (r_pos >= (q_pos-align_length)) & (r_pos <= (q_pos+args.L)):
                        if not nr in output_dict.keys():
                            if q_pos == r_pos:
                                cov_bases = align_length
                            else:
                                cov_bases = align_length - abs(r_pos-q_pos)
                            output_dict[nr] = (1,cov_bases)
                        else:
                            print("[Warning] Multiple mappings per query")
                    else:
                        # print("[Fwd]", line[0:4])
                        counter_fw += 1
                        if not nr in output_dict.keys():
                            output_dict[nr] = (0,0)
                        else:
                            print("[Warning] Multiple mappings per query")

                elif line[1] == "16":  # primary alignment rev
                    align_length = len(line[9])
                    nr, chr, dir, q_pos, mut = line[0].split("|")
                    # print(nr, chr, dir, q_pos, mut)
                    q_pos = int(q_pos)
                    r_pos = int(line[3])
                    if (dir=="rev") & (line[2] == chr) & (r_pos >= (q_pos-align_length)) & (r_pos <= (q_pos+args.L)):
                        if not nr in output_dict.keys():
                            if q_pos == r_pos:
                                cov_bases = align_length
                            else:
                                cov_bases = align_length - abs(r_pos-q_pos)
                            output_dict[nr] = (1,cov_bases)
                        else:
                            print("[Warning] Multiple mappings per query")
                    else:
                        # print("[Rev]", line[0:4])
                        counter_rev += 1
                        if not nr in output_dict.keys():
                            output_dict[nr] = (0,0)
                        else:
                            print("[Warning] Multiple mappings per query")

    # print(output_dict)
    mapped, covered = zip(*output_dict.values())
    print(sum(mapped)/len(output_dict))
    print(sum(covered)/args.N)
    print(counter_fw, counter_rev, len(output_dict))




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_file', type=str, help='Name of SAM output file from minimap2')
    parser.add_argument('--L', type=int, default=10000, help='Length of each query read')
    parser.add_argument('--N', type=int, default=100000, help='Number of query reads')
    args = parser.parse_args()

    main(args)
