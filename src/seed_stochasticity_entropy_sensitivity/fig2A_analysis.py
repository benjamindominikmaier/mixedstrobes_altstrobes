import os,sys
import argparse

import random
import math
from collections import defaultdict

# import numpy as np
from math import log, e

import indexing_Maier_altstrobes_new_multi as indexing
# try:
#     import matplotlib
#     matplotlib.use('Agg')
#     import matplotlib.pyplot as plt
# except (ImportError, RuntimeError):
#     print("COULD not import matplotlib")

# import seaborn as sns
# import pandas as pd
# from matplotlib import pyplot

# def entropy(labels, base=None):
#   """ Computes entropy of label distribution. """

#   n_labels = len(labels)

#   if n_labels <= 1:
#     return 0

#   value,counts = np.unique(labels, return_counts=True)
#   probs = counts / n_labels
#   # print(probs)
#   n_classes = np.count_nonzero(probs)

#   if n_classes <= 1:
#     return 0

#   ent = 0.

#   # Compute entropy
#   base = e if base is None else base
#   for i in probs:
#     ent -= i * log(i, base)

#   return ent

def entropy(labels, base=None):
    return 0

def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def plot_hist(nr_matches, plot_file_prefix):
    plt.hist(nr_matches, density=False) #, bins=100, range=[0, 100])  # density=False would make counts
    # plt.yscale('log', nonposy='clip')
    plt.ylabel('Probability')
    plt.xlabel('Number of matches')
    plt.savefig(plot_file_prefix +'.pdf' ,bbox_inches='tight')
    plt.clf()

def compute_frac_matches_per_mut_rate(args, N_sim, outfile):
    k = args.k
    w_min = args.w_min
    w_max = args.w_max
    w = w_max + k//2 - 1
    L = 4*w # for W=50, guarnatee 50 seeds sampled and check P(N_m = 0) within [0,W] window for 2*m mutations in the string pair
    # pre generate list of sting pairs for controlled comparison
    seq_pairs = defaultdict(list)
    freq = [0.02, 0.05, 0.08, 0.1, 0.15]
    #N_m = [int(2*w_max*freq[0]),int(2*w_max*freq[1]),int(2*w_max*freq[2]), int(2*w_max*freq[3]), int(2*w_max*freq[4])] # [2, 5, 10, 15]  #
    #N_m = [i for i in range(1,21)]
    N_m = [i for i in range(1,26)]
    # m_to_freq = { N_m[i] : freq[i] for i in range(len(N_m)) }
    print(N_m)
    for m in N_m:
        seq_pairs[m] = []
        for n in range(N_sim):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])
            muts = set(random.sample(range(2*w),m)) # only sample mutations within the 2W region to control mutation rate in the region according to theory
            subs = random.uniform(0, 1)
            seq2 = "".join([seq1[i] if i not in muts else random.choices([reverse_complement(seq1[i]), '',  seq1[i] + random.choice("ACGT")],  weights = [subs, subs/2,subs/2] )[0] for i in range(len(seq1))])
            seq_pairs[m].append((seq1,seq2))
        print(muts)
        print(len(seq_pairs[m]))


    # k-mers

    # calc Nr matches
    cumulative_frac = []
    for m in N_m:
        # print('for m=', m)
        had_a_match = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[m][n]
            kmers_pos1 = indexing.kmers(seq1, k, 1)
            kmers_pos2 = indexing.kmers(seq2, k, 1)

            # check the seeds within the window
            seeds1 = [s for i, s in kmers_pos1.items() if i <= w and i % args.thinning == 0]
            seeds2 = [s for i, s in kmers_pos2.items() if i <= w and i % args.thinning == 0]
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)
            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)

        entr = entropy(had_a_match,base=2)
        frac_matches = sum(had_a_match)/N_sim
        stddev_frac_matches = math.sqrt(sum([(n - frac_matches)**2 for n in had_a_match])/N_sim)
        cumulative_frac.append(frac_matches)

        print("k-mers: frac_matches: {0}, stddev_frac_matches: {1}, err rate:{2}, Est entropy:{3}".format(frac_matches, stddev_frac_matches, m, entr) )
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(w_min, w_max, 'k-mers', frac_matches, stddev_frac_matches, m))
        # plot_hist(had_a_match, args.outplot_prefix+'_kmers_{0}'.format(m) )
    print('k-mers sum frac:', sum(cumulative_frac))
    print()

    # Minstrobes
    cumulative_frac = []
    for m in N_m:
        # print(s, 'for m=', m)
        had_a_match = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[m][n]
            minstrobes_pos1 = indexing.minstrobes(seq1, k, w_min, w_max, 1)
            # print(minstrobes_pos1)
            minstrobes_pos2 = indexing.minstrobes(seq2, k, w_min, w_max, 1)

            # check the seeds within the window
            seeds1 = [s for i_tuple, s in minstrobes_pos1.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            seeds2 = [s for i_tuple, s in minstrobes_pos2.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            # print(len(seeds1))
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)
            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)

        entr = entropy(had_a_match,base=2)
        frac_matches = sum(had_a_match)/N_sim
        stddev_frac_matches = math.sqrt(sum([(n - frac_matches)**2 for n in had_a_match])/N_sim)
        cumulative_frac.append(frac_matches)

        print("minstrobes: frac_matches: {0}, stddev_frac_matches: {1}, err rate:{2}, Est entropy:{3}".format(frac_matches, stddev_frac_matches, m, entr) )
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(w_min, w_max, 'minstrobes', frac_matches, stddev_frac_matches, m))
        # plot_hist(had_a_match, args.outplot_prefix+'_randstrobes_{0}'.format(m) )
    print('Minstrobes sum frac:', sum(cumulative_frac))
    print()

    # Hybridstrobes
    cumulative_frac = []
    for m in N_m:
        # print(s, 'for m=', m)
        had_a_match = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[m][n]
            hybridstrobes_pos1 = indexing.hybridstrobes(seq1, k, w_min, w_max, 1)
            # print(hybridstrobes_pos1)
            hybridstrobes_pos2 = indexing.hybridstrobes(seq2, k, w_min, w_max, 1)

            # check the seeds within the window
            seeds1 = [s for i_tuple, s in hybridstrobes_pos1.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            seeds2 = [s for i_tuple, s in hybridstrobes_pos2.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            # print(len(seeds1))
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)
            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)

        entr = entropy(had_a_match,base=2)
        frac_matches = sum(had_a_match)/N_sim
        stddev_frac_matches = math.sqrt(sum([(n - frac_matches)**2 for n in had_a_match])/N_sim)
        cumulative_frac.append(frac_matches)

        print("hybridstrobes: frac_matches: {0}, stddev_frac_matches: {1}, err rate:{2}, Est entropy:{3}".format(frac_matches, stddev_frac_matches, m, entr) )
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(w_min, w_max, 'hybridstrobes', frac_matches, stddev_frac_matches, m))
        # plot_hist(had_a_match, args.outplot_prefix+'_randstrobes_{0}'.format(m) )
    print('Hybridstrobes sum frac:', sum(cumulative_frac))
    print()


    # Randstrobes
    cumulative_frac = []
    k_s = k//2
    k_l = k//2

    # calc P match
    for m in N_m:
        # print(s, 'for m=', m)
        had_a_match = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[m][n]
            altstrobes_pos1 = indexing.altstrobes_generalized(seq1, k_s, k_l, w_min, w_max, 1)
            # print(altstrobes_pos1)
            altstrobes_pos2 = indexing.altstrobes_generalized(seq2, k_s, k_l, w_min, w_max, 1)

            # check the seeds within the window
            seeds1 = [s for i_tuple, s in altstrobes_pos1.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            seeds2 = [s for i_tuple, s in altstrobes_pos2.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            # print(len(seeds1))
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)
            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)

        entr = entropy(had_a_match,base=2)
        frac_matches = sum(had_a_match)/N_sim
        stddev_frac_matches = math.sqrt(sum([(n - frac_matches)**2 for n in had_a_match])/N_sim)
        cumulative_frac.append(frac_matches)

        print("randstrobes: frac_matches: {0}, stddev_frac_matches: {1}, err rate:{2}, Est entropy:{3}".format(frac_matches, stddev_frac_matches, m, entr) )
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(w_min, w_max, 'randstrobes', frac_matches, stddev_frac_matches, m))
        # plot_hist(had_a_match, args.outplot_prefix+'_randstrobes_{0}'.format(m) )
    print('Randstrobes sum frac:', sum(cumulative_frac))
    print()

    # Altstrobes 10/20
    k_s = 10
    k_l = 20
    cumulative_frac = []
    # calc P match
    for m in N_m:
        # print(s, 'for m=', m)
        had_a_match = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[m][n]
            altstrobes_pos1 = indexing.altstrobes_generalized(seq1, k_s, k_l, w_min, w_max, 1)
            # print(altstrobes_pos1)
            altstrobes_pos2 = indexing.altstrobes_generalized(seq2, k_s, k_l, w_min, w_max, 1)

            # check first W seeds
            seeds1 = [s for i_tuple, s in altstrobes_pos1.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            seeds2 = [s for i_tuple, s in altstrobes_pos2.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)

            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)

        entr = entropy(had_a_match,base=2)
        frac_matches = sum(had_a_match)/N_sim
        stddev_frac_matches = math.sqrt(sum([(n - frac_matches)**2 for n in had_a_match])/N_sim)
        cumulative_frac.append(frac_matches)

        print("Altstrobes: frac_matches: {0}, stddev_frac_matches: {1}, err rate:{2}, Est entropy:{3}".format(frac_matches, stddev_frac_matches, m, entr) )
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(w_min, w_max, 'altstrobes', frac_matches, stddev_frac_matches, m))
        # plot_hist(had_a_match, args.outplot_prefix+'_altstrobes_{0}'.format(m) )
    print('Altstrobes sum frac:', sum(cumulative_frac))
    print()

    # Mixedstrobes 80-20
    p_strobe = 8/10
    k_s = k//2
    k_l = k//2
    cumulative_frac = []
    # calc P match
    for m in N_m:
        # print(s, 'for m=', m)
        had_a_match = []
        all_var = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[m][n]
            mixedstrobes_pos1 = indexing.mixedrandstrobes(seq1, k, w_min, w_max, 1, strobe_fraction=p_strobe)
            # print(mixedstrobes_pos1)
            mixedstrobes_pos2 = indexing.mixedrandstrobes(seq2, k, w_min, w_max, 1, strobe_fraction=p_strobe)

            # check first W seeds
            seeds1 = [s for i_tuple, s in mixedstrobes_pos1.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            seeds2 = [s for i_tuple, s in mixedstrobes_pos2.items() if i_tuple[0] <= w and i_tuple[0] % args.thinning == 0]
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)

            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)

        entr = entropy(had_a_match,base=2)
        frac_matches = sum(had_a_match)/N_sim
        stddev_frac_matches = math.sqrt(sum([(n - frac_matches)**2 for n in had_a_match])/N_sim)
        cumulative_frac.append(frac_matches)

        print("Mixedstrobes: frac_matches: {0}, stddev_frac_matches: {1}, err rate:{2}, Est entropy:{3}".format(frac_matches, stddev_frac_matches, m, entr) )
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(w_min, w_max, 'mixedstrobes', frac_matches, stddev_frac_matches, m))
        # plot_hist(had_a_match, args.outplot_prefix+'_mixedstrobes_{0}'.format(m) )
    print('Mixedstrobes sum frac:', sum(cumulative_frac))
    print()


    # Altstrobes-multisized k_lower = 4
    k_s = 4
    cumulative_frac = []
    # calc P match
    for m in N_m:
        # print(s, 'for m=', m)
        had_a_match = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[m][n]
            altstrobes_pos1 = indexing.altstrobes_continous(seq1, k, w_min, w_max, 1, k_boundary=k_s)
            # print([i_tuple for i_tuple, s in altstrobes_pos1.items() if i_tuple[0][0] <= w] )
            # print([i_tuple[0][0] for i_tuple, s in altstrobes_pos1.items() if i_tuple[0][0] <= w] )
            altstrobes_pos2 = indexing.altstrobes_continous(seq2, k, w_min, w_max, 1, k_boundary=k_s)

            # check first W seeds
            seeds1 = [s for i_tuple, s in altstrobes_pos1.items() if i_tuple[0][0] <= w and i_tuple[0][0] % args.thinning == 0 ]
            seeds2 = [s for i_tuple, s in altstrobes_pos2.items() if i_tuple[0][0] <= w and i_tuple[0][0] % args.thinning == 0 ]
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)

            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)

        entr = entropy(had_a_match,base=2)
        frac_matches = sum(had_a_match)/N_sim
        stddev_frac_matches = math.sqrt(sum([(n - frac_matches)**2 for n in had_a_match])/N_sim)
        cumulative_frac.append(frac_matches)

        print("Altstrobes-multi: frac_matches: {0}, stddev_frac_matches: {1}, err rate:{2}, Est entropy:{3}".format(frac_matches, stddev_frac_matches, m, entr) )
        outfile.write("{0},{1},{2},{3},{4},{5}\n".format(w_min, w_max, 'multistrobes', frac_matches, stddev_frac_matches, m))
        # plot_hist(had_a_match, args.outplot_prefix+'_altstrobes_{0}'.format(m) )
    print('Altstrobes sum frac:', sum(cumulative_frac))
    print()

    return

def plot_p_match(args, input_csv, plot_file_prefix , x="k_2/k"):

    palette = {
    'k-mers': 'tab:blue',
    'minstrobes': 'tab:grey',
    'hybridstrobes': 'gold',
    'randstrobes': 'tab:red',
    'mixedstrobes' : 'tab:purple',
    'altstrobes' : 'tab:green',
    'multistrobes' : 'black'
    # "bwa_mem2" : 'pink',
    # "strobealign_master" : 'cyan',
    # "strobealign_master_preindexed" : 'gold'
    }
    methods = ['k-mers', 'minstrobes', 'hybridstrobes', 'randstrobes', 'mixedstrobes', 'altstrobes', 'multistrobes'] #, "strobealign_mixed"] "urmap",

    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    indata = pd.read_csv(input_csv)
    # sns.set_palette("pastel", 8, .75)
    print(sns.color_palette())
    sns.set_palette("Set1")
    g = sns.relplot(data=indata, x="m", y="N_m", hue="seeding", kind='line', facet_kws=dict(sharex=False),
                    hue_order = methods,
                    palette=palette) #s=80,
    # g.set(yscale="log")
    # g.set(ylim=(0.0001, 1.2))
    plt.xticks([1, 5, 10,15, 20])
    g.set_ylabels("P(N_m > 0)")
    plt.savefig(plot_file_prefix +'_frac_matches.pdf' ,bbox_inches='tight')
    plt.close()

def main(args):
    outfile = open(args.outcsv, 'w')
    outfile.write("w_min,w_max,seeding,N_m,stddev_est,m\n") # N is nr matches, N_mut is nr muations
    print('Using N_SIM:', args.N_SIM)
    compute_frac_matches_per_mut_rate(args, args.N_SIM, outfile)
    outfile.close()
    plot_p_match(args, outfile.name, args.outplot_prefix)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    parser.add_argument('k', type=int, default=15, help='k-mer size')
    parser.add_argument('w_min', type=int, default=25, help='Window size')
    parser.add_argument('w_max', type=int, default=50, help='Window size')
    # parser.add_argument('f', type=str, default="2/3", help='relative size between strobes within strobemer')
    parser.add_argument('outcsv', type=str,  default=None, help='output CSV file for plotting.')
    parser.add_argument('outplot_prefix', type=str,  default=None, help='output pdf file prefix.')
    parser.add_argument('--N_SIM', type=int, default=1000, help='Number of simulations for p_matches experiment')
    parser.add_argument('--thinning', type=int, default=1, help='Thinning level, 1/X. Default X=1 no thinning. ')

    # parser.set_defaults(which='main')
    args = parser.parse_args()



    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit()

    main(args)
