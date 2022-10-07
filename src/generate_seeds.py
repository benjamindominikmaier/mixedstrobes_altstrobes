#!/usr/bin/env python3.9
# -*- coding: utf-8 -*
import os
import sys
import argparse
import random
from itertools import chain
import more_itertools as mit


def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break

def mutate_sequences(args, seq1, mut_freq):
    """
    """
    # controlled or random experiment
    if args.experiment_type == 'only_subs':
        muts = set(random.sample(range(len(seq1)), int(args.L*mut_freq)))
        seq2 = "".join([
            seq1[i] if i not in muts
            else random.choice([reverse_complement(seq1[i])])
            for i in range(len(seq1))
        ])
    elif args.experiment_type == 'controlled':
        # muts = set(range(15,L,15)) # every 15th nt for figure 2 only!
        muts = set(range(20, args.L, 20))
        seq2 = "".join([
            seq1[i] if i not in muts
            else random.choice(['', reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")])
            for i in range(len(seq1))
        ])
    elif args.experiment_type == 'all':
        muts = set(random.sample(range(len(seq1)), int(args.L*mut_freq)))
        seq2 = "".join([
            seq1[i] if i not in muts
            else random.choice(['', reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")])
            for i in range(len(seq1))
        ])
    else:
        print("Wrong experiment label specified")
        sys.exit()
    return seq2


def print_sim_reads_to_file(args, reads, mut_freq) -> None:
    """
    Print matches to fasta-file.

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
    if mut_freq > 0:
        outfile = open(os.path.join(args.outfolder, args.prefix + "_" + str(mut_freq) + "_mutated.txt"), 'w')
    else:
        outfile = open(os.path.join(args.outfolder, args.prefix), 'w')

    for nr, read in enumerate(reads, 1):
        name, dir, pos, seq = read
        outfile.write(">{0}|{1}|{2}|{3}|{4}\n{5}\n".format(nr, name, dir, pos, mut_freq, seq))


def reverse_complement(string):
    """
    """
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def main(args):
    """
    """
    fasta_records = {}
    len_fasta_records = {}
    print("Start processing fasta records from {0}".format(args.reference))
    for nr, record in enumerate(readfq(open(args.reference, 'r')), 1):
        name, fasta_record = record
        fasta_record = fasta_record[0].upper()
        # print("Process fasta_entry:", name, len(fasta_record))
        if len(fasta_record) >= args.L:
            if args.tolerate_N:
                len_fasta_records[name] = len(fasta_record)-args.L
            else:
                len_fasta_records[name] = len(fasta_record.strip("N"))-args.L
    reference_seed_distribution = random.choices(list(len_fasta_records.keys()), weights = list(len_fasta_records.values()), k = args.nr_exp)

    print("Simulate {0} fasta reads from {1}".format(args.nr_exp, args.reference))
    sim_reads = []
    for nr, record in enumerate(readfq(open(args.reference, 'r')), 1):
        name, fasta_record = record
        fasta_record = fasta_record[0].upper()
        n = reference_seed_distribution.count(name)
        if n > 0:
            if args.tolerate_N:
                read_starts = random.sample(range(len(fasta_record)-args.L), n)
                for read_start in read_starts:
                    seq = fasta_record[read_start: read_start+args.L]
                    if args.rev_comp & random.choice([0,1]):
                        sim_reads.append((name, "rev", read_start+1, reverse_complement(seq)))
                    else:
                        sim_reads.append((name, "forw", read_start+1, seq))
                    if len(sim_reads) % 1000 == 0:
                        print("Simulated {0} fasta reads".format(len(sim_reads)))
            else:
                fasta_record_stripped = fasta_record.lstrip("N")
                shift = len(fasta_record)-len(fasta_record_stripped)
                # print("Removed {0} leading Ns".format(shift))
                fasta_record_stripped = fasta_record.rstrip("N")
                reads_seeded = 0
                abort = 0
                while reads_seeded < n:
                    read_start = random.sample(range(len(fasta_record_stripped)-args.L), 1)[0]
                    seq = fasta_record_stripped[read_start: read_start+args.L]
                    if not "N" in seq:
                        if args.rev_comp & random.choice([0,1]):
                            sim_reads.append((name, "rev", read_start+1, reverse_complement(seq)))
                        else:
                            sim_reads.append((name, "forw", read_start+1, seq))
                        reads_seeded += 1
                        if len(sim_reads) % 1000 == 0:
                            print("Simulated {0} fasta reads".format(len(sim_reads)))
                    else:
                        abort += 1
                        if abort > 100*n:
                            print("\t[Warning] Unable to seed enough reads from fasta record {0} due to too many unknown nucleotides. {1}/{2} reads were seeded from sequence".format(name, reads_seeded, n))
                            break

    print("\tWrite simulated fasta reads to file.")
    print_sim_reads_to_file(args, sim_reads, 0)
    for mut_freq in args.mut_freqs:
        print("Mutate simulated fasta reads (mutuation rate = {0})".format(mut_freq))
        mut_sim_reads = []
        for sim_read in sim_reads:
            name, dir, pos, seq = sim_read
            mut_sim_reads.append((name, dir, pos, mutate_sequences(args, seq, mut_freq)))
        print("\tWrite mutated fasta reads to file.")
        print_sim_reads_to_file(args, mut_sim_reads, mut_freq)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--L', type=int, default=10000, help='Length of simulated sequences')
    parser.add_argument('--nr_exp', type=int, default=100000, help='Number of simulated experiments')
    parser.add_argument('--experiment_type', type=str, default="all", help='experiment type choose between "all", "controlled or "only_subs"')
    parser.add_argument('--mut_freqs', type=list, default=[0.01, 0.05, 0.1], help='mutation frequencies [0.01, 0.05, 0.1]')
    parser.add_argument('--reference', type=str, help='path to fasta file with reference genome to generate reads from')
    parser.add_argument('--rev_comp', action='store_true', help='Sampling reads from both strands')
    parser.add_argument('--outfolder', type=str,  default="../output/simulated_reads", help='Folder to output FA read file.')
    parser.add_argument('--prefix', type=str,  default="reads", help='Filename prefix (default "matches").')
    parser.add_argument('--tolerate_N', action='store_true', help="Allow Ns in simulated reads")
    args = parser.parse_args()

    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
