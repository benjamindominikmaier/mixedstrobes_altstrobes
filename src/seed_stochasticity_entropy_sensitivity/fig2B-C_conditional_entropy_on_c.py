import os,sys
import argparse

import random
import math
from collections import defaultdict

import indexing_Maier_altstrobes as indexing

# try:
#     import matplotlib
#     matplotlib.use('Agg')
#     import matplotlib.pyplot as plt
# except (ImportError, RuntimeError):
#     print("COULD not import matplotlib")

# import seaborn as sns
# import pandas as pd
# from matplotlib import pyplot


"""
Main theorem:
Assume a positon j in S and a seed construct function f(\cdot, k, W, S).
Let p(s_i) be the probability that seed s_i produced by f covers position j.

The most efficient seed sampling protocol f will
be the one that maximizes entropy
H(X) = - \sum_{i=0}^{|S|- W +1} p(s_i) log p(s_i)
= - \sum_{i=j-W}^{j} p(s_i) log p(s_i)
"""

def plot_entropy(args, input_csv, plot_file_prefix , x="k_2/k"):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    indata = pd.read_csv(input_csv)
    indata_entropy = indata.loc[indata['ymetric'].str.contains('Entropy')]
    print(indata_entropy)

    g = sns.relplot(data=indata_entropy, x="x", y="y", hue="type", style="w_min", kind="line", col="analysis", facet_kws=dict(sharex=False)) #s=80,
    g.axes[0,0].set_xlim(0,31)
    g.axes[0,0].set_xlabel('Length k_s')
    g.axes[0,1].set_xlim(-0.1,1.1)
    g.axes[0,1].set_xlabel('fraction randstrobes')
    #g.axes[0,0].set_ylim(-1,23)
    axes = g.axes
    # g.set(yscale="log")
    # g.set( yticks= [9] + [i for i in range(10,99,10)] + [i for i in range(100,999,100)])
    # print(g._legend.texts)
    # g._legend.texts[0].set_text(r'$k$-mers')

    # print(axes)
    axes = g.axes.flatten()
    # axes[0].set_title("With 1 nbr seed")
    # axes[1].set_title("With 5 nbr seeds")
    # axes[2].set_title("With {0} nbr seeds".format(args.w_max + args.k) )
    axes[0].set_ylabel("Entropy")
    # /P(N_m)
    #axes[3].set_ylabel("% Expected overlapping nt")

    plt.savefig(plot_file_prefix +'_entropy.pdf' ,bbox_inches='tight')
    plt.close()


def plot_p_match(args, input_csv, plot_file_prefix , x="k_2/k"):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    indata = pd.read_csv(input_csv)
    indata_p_match = indata.loc[indata['ymetric'].str.contains('p_match')]
    print(indata_p_match)

    g = sns.relplot(data=indata_p_match, x="x", y="y", hue="type", style="w_min", kind="line", col="analysis", facet_kws=dict(sharex=False))  #s=80,
    g.axes[0,0].set_xlim(0,31)
    g.axes[0,0].set_xlabel('Length k_s')
    g.axes[0,1].set_xlim(-0.1,1.1)
    g.axes[0,1].set_xlabel('fraction randstrobes')
    #g.axes[0,0].set_ylim(-0.1,1.0)
    axes = g.axes
    # g.set(yscale="log")
    # g.set( yticks= [9] + [i for i in range(10,99,10)] + [i for i in range(100,999,100)])
    # print(g._legend.texts)
    # g._legend.texts[0].set_text(r'$k$-mers')

    # print(axes)
    axes = g.axes.flatten()
    # axes[0].set_title("With 1 nbr seed")
    # axes[1].set_title("With 5 nbr seeds")
    # axes[2].set_title("With {0} nbr seeds".format(args.w_max + args.k) )
    axes[0].set_ylabel("P(N_m)")
    # /P(N_m)
    #axes[3].set_ylabel("% Expected overlapping nt")

    plt.savefig(plot_file_prefix +'_p_match.pdf' ,bbox_inches='tight')
    plt.close()


def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def compute_p_match(args, N_sim, outfile):
    k = args.k
    w_min = args.w_min
    w_max = args.w_max
    W = w_max + k//2 - 1
    L = 4*W # for W=64, guarnatee 64 seeds sampled and check P(N_m = 0) within [0,W] window for 2*m mutations in the string pair
    # pre generate list of sting pairs for controlled comparison
    seq_pairs = defaultdict(list)
    # freq = [0.02, 0.05, 0.08, 0.1, 0.15]
    # N_m = [int(2*w_max*freq[0]),int(2*w_max*freq[1]),int(2*w_max*freq[2]), int(2*w_max*freq[3]), int(2*w_max*freq[4])] # [2, 5, 10, 15]  #
    # m_to_freq = { N_m[i] : freq[i] for i in range(len(N_m)) }
    # N_m = [5, 10, 15]

    MAX_MUT = int(args.max_mut_freq*2*W)
    N_m = [i for i in range(1,MAX_MUT+1)]

    # N_m = [i for i in range(1,61)]
    # N_m = [i for i in range(26,31)]
    print(N_m, W)

    ######## initialization  ##########

    kmer_p_matches = []

    mixedstrobes_p_matches = {}
    for f in range(0, 11):
        frac_strobes = f/10
        mixedstrobes_p_matches[frac_strobes] = []

    altstrobe_p_matches = {}
    for k_s in range(1, k//2+1):
        altstrobe_p_matches[k_s] = []

    multistrobe_p_matches = {}
    for k_s in range(1, k//2+1):
        multistrobe_p_matches[k_s] = []
    ###################################


    for m in N_m:
        print('m=',m)
        seq_pairs = []
        for n in range(N_sim):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])
            muts = set(random.sample(range(2*W),m)) # only sample mutations within the 2W region to control mutation rate in the region according to theory
            subs = random.uniform(0, 1)
            seq2 = "".join([seq1[i] if i not in muts else random.choices([reverse_complement(seq1[i]), '',  seq1[i] + random.choice("ACGT")],  weights = [subs, subs/2,subs/2] )[0] for i in range(len(seq1))])
            seq_pairs.append((seq1,seq2))
        print(muts)
        print(len(seq_pairs))


        had_a_match = []
        for n in range(N_sim):
            # print('for n=', n)
            seq1 ,seq2 = seq_pairs[n]
            kmers_pos1 = indexing.kmers(seq1, k, 1)
            kmers_pos2 = indexing.kmers(seq2, k, 1)

            # check first W seeds
            seeds1 = [s for i, s in kmers_pos1.items() if i < W]
            seeds2 = [s for i, s in kmers_pos2.items() if i < W]
            matches = set(seeds1) & set(seeds2)
            n_hits = len(matches)
            if n_hits > 0:
                had_a_match.append(1)
            else:
                had_a_match.append(0)


        p_match = sum(had_a_match)/N_sim
        kmer_p_matches.append(p_match)
        # mean_error_p_match = math.sqrt(sum([(x - p_match)**2 for x in had_a_match]) / N_sim) / math.sqrt(N_sim)
        print("k-mers: p_match: {0}, m:{1}".format(p_match, m) )


        # Altstrobes
        for s in range(1, k//2+1):
            k_s = s
            k_l = k-s

            # print(s, 'for m=', m)
            had_a_match = []
            all_var = []
            for n in range(N_sim):
                # print('for n=', n)
                seq1 ,seq2 = seq_pairs[n]
                altstrobes_pos1 = indexing.altstrobes_generalized(seq1, k_s, k_l, w_min, w_max, 1)
                # print(altstrobes_pos1)
                altstrobes_pos2 = indexing.altstrobes_generalized(seq2, k_s, k_l, w_min, w_max, 1)

                # check first W seeds
                seeds1 = [s for i_tuple, s in altstrobes_pos1.items() if i_tuple[0] < W]
                seeds2 = [s for i_tuple, s in altstrobes_pos2.items() if i_tuple[0] < W]
                matches = set(seeds1) & set(seeds2)
                n_hits = len(matches)

                if n_hits > 0:
                    had_a_match.append(1)
                else:
                    had_a_match.append(0)

            p_match = sum(had_a_match)/N_sim
            # mean_error_p_match = math.sqrt(sum([(x - p_match)**2 for x in had_a_match]) / N_sim) / math.sqrt(N_sim)
            altstrobe_p_matches[k_s].append(p_match)
            print("altstrobes-{0}: p_match: {1}, m:{2}".format(k_s, p_match, m) )


        # Mixedstrobes
        # Assumes  k_l = k_s  < w_min (the strobes are smaller than w_min)
        for f in range(0, 11):
            p_strobe = f/10
            k_s = k//2
            k_l = k//2
            # print(s, 'for m=', m)
            had_a_match = []
            all_var = []
            for n in range(N_sim):
                # print('for n=', n)
                seq1 ,seq2 = seq_pairs[n]
                mixedstrobes_pos1 = indexing.mixedrandstrobes(seq1, k, w_min, w_max, 1, strobe_fraction=p_strobe)
                # kmers1 = len([s for i_tuple, s in mixedstrobes_pos1.items() if (i_tuple[0] < W) and (i_tuple[0]+k_s == i_tuple[1])])
                mixedstrobes_pos2 = indexing.mixedrandstrobes(seq2, k, w_min, w_max, 1, strobe_fraction=p_strobe)
                # kmers2 = len([s for i_tuple, s in mixedstrobes_pos2.items() if (i_tuple[0] < W) and (i_tuple[0]+k_s == i_tuple[1])])
                # print('M {0}: k-mers1: {1}, k-mers2: {2}'.format(p_strobe, kmers1, kmers2))

                # check first W seeds
                seeds1 = [s for i_tuple, s in mixedstrobes_pos1.items() if i_tuple[0] < W]
                seeds2 = [s for i_tuple, s in mixedstrobes_pos2.items() if i_tuple[0] < W]
                matches = set(seeds1) & set(seeds2)
                n_hits = len(matches)

                if n_hits > 0:
                    had_a_match.append(1)
                else:
                    had_a_match.append(0)

            p_match = sum(had_a_match)/N_sim
            # mean_error_p_match = math.sqrt(sum([(x - p_match)**2 for x in had_a_match]) / N_sim) / math.sqrt(N_sim)
            mixedstrobes_p_matches[p_strobe].append(p_match)
            print("mixedstrobes-{0}: p_match: {1}, m:{2}".format(p_strobe, p_match, m) )

        print()


        # Multistrobes
        for k_s in range(1,k//2+1):
            had_a_match = []
            for n in range(N_sim):
                # print('for n=', n)
                seq1 ,seq2 = seq_pairs[n]
                altstrobes_pos1 = indexing.altstrobes_continous(seq1, k, w_min, w_max, 1, k_boundary=k_s)
                # print([i_tuple for i_tuple, s in altstrobes_pos1.items() if i_tuple[0][0] <= w] )
                # print([i_tuple[0][0] for i_tuple, s in altstrobes_pos1.items() if i_tuple[0][0] <= w] )
                altstrobes_pos2 = indexing.altstrobes_continous(seq2, k, w_min, w_max, 1, k_boundary=k_s)

                # check first W seeds
                seeds1 = [s for i_tuple, s in altstrobes_pos1.items() if i_tuple[0][0] < W ]
                seeds2 = [s for i_tuple, s in altstrobes_pos2.items() if i_tuple[0][0] < W ]
                # print(seeds1)
                # print(k_s, len(seeds1), len(seeds2))
                matches = set(seeds1) & set(seeds2)
                n_hits = len(matches)

                if n_hits > 0:
                    had_a_match.append(1)
                else:
                    had_a_match.append(0)

            p_match = sum(had_a_match)/N_sim
            # mean_error_p_match = math.sqrt(sum([(x - p_match)**2 for x in had_a_match]) / N_sim) / math.sqrt(N_sim)
            multistrobe_p_matches[k_s].append(p_match)
            print("multistrobes-{0}: p_match: {1}, m:{2}".format(k_s, p_match, m) )

        print()

    ###################### OUTPUT ########################

    print('FINAL')

    print("k-mers: p_match: {0}, m:{1}, mean error est: {2}".format(sum(kmer_p_matches), m, '-') )
    outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'k-mers', k, m, sum(kmer_p_matches), '-', 'altstrobes-analysis', 'p_match'))


    for s in range(1, k//2+1):
        k_s = s
        k_l = k-s
        print("Altstrobes: ({0},{1}), p_match: {2}, m:{3}, mean error est: {4}".format( k_s, k_l, sum(altstrobe_p_matches[k_s]), m, '-' ))
        if k_s == k_l:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', k_s, m, sum(altstrobe_p_matches[k_s]), '-', 'altstrobes-analysis', 'p_match'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'altstrobes', k_s, m, sum(altstrobe_p_matches[k_s]), '-', 'altstrobes-analysis', 'p_match'))


    for f in range(0, 11):
        p_strobe = f/10
        k_s = k//2
        k_l = k//2
        print("mixedstrobes-{0}% p_match: {1}, m:{2}, mean error est: {3}".format(round(100*p_strobe,0), sum(mixedstrobes_p_matches[p_strobe]), m, '-' ) )
        if p_strobe == 0:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'k-mers', p_strobe,  m, sum(mixedstrobes_p_matches[p_strobe]), '-', 'mixedstrobes-analysis', 'p_match'))
        elif p_strobe == 1.0:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', p_strobe,  m, sum(mixedstrobes_p_matches[p_strobe]), '-', 'mixedstrobes-analysis', 'p_match'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'mixed-randstrobes', p_strobe,  m, sum(mixedstrobes_p_matches[p_strobe]), '-', 'mixedstrobes-analysis', 'p_match'))

    for s in range(1, k//2+1):
        k_s = s
        k_l = k-s
        print("Multistrobes: ({0},{1}), p_match: {2}, m:{3}, mean error est: {4}".format( k_s, k_l, sum(multistrobe_p_matches[k_s]), m, '-' ))
        if k_s == k_l:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', k_s, m, sum(multistrobe_p_matches[k_s]), '-', 'multistrobes-analysis', 'p_match'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'multistrobes', k_s, m, sum(multistrobe_p_matches[k_s]), '-', 'multistrobes-analysis', 'p_match'))


    print()

    ######################################################



    # # Randstrobes order 3
    # p_matches = []
    # # k_strobe = k//3
    # w_min_r3 = k//3 + 1
    # w_max_r3 = w_max // 2 + 10
    # for m in N_m:
    #     # print('for m=', m)
    #     had_a_match = []
    #     for n in range(N_sim):
    #         # print('for n=', n)
    #         seq1 ,seq2 = seq_pairs[n]
    #         randstrobes_pos1 = indexing.randstrobes(seq1, k, w_min_r3, w_max_r3, 1, order=3)
    #         # print(altstrobes_pos1)
    #         randstrobes_pos2 = indexing.randstrobes(seq2, k, w_min_r3, w_max_r3, 1, order=3)
    #         #print([i_tuple for i_tuple, s in randstrobes_pos1.items() if i_tuple[0] < W] )
    #         # check first W seeds
    #         seeds1 = [s for i_tuple, s in randstrobes_pos1.items() if i_tuple[0] < W]
    #         seeds2 = [s for i_tuple, s in randstrobes_pos2.items() if i_tuple[0] < W]
    #         # print(len(seeds1))
    #         matches = set(seeds1) & set(seeds2)
    #         n_hits = len(matches)
    #         if n_hits > 0:
    #             had_a_match.append(1)
    #         else:
    #             had_a_match.append(0)


    #     p_match = sum(had_a_match)/N_sim
    #     print("Randstrobes Order3: p_match: {0}, m:{1}, mean error est: {2}".format(p_match, m, '-') )

    #     p_matches.append(p_match)
    #     # mean_error_p_match = math.sqrt(sum([(x - p_match)**2 for x in had_a_match]) / N_sim) / math.sqrt(N_sim)
    # print("Randstrobes Order3: p_match: {0}, m:{1}, mean error est: {2}".format(sum(p_matches), m, '-') )
    # outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'k-mers', k, m, sum(p_matches), '-', 'altstrobes-analysis', 'p_match'))

    return




def P_X_i_given_Y_j(i, j, s1, s2, w_min, w_max):
    """
        s1: strobe 1 length.
        s2: strobe 2 length.

        j is in [1,k]
        This function assumes w_min > s2 for the probabilities to be correct
    """

    # ################################################################################
    # #### Below if-else takes care of impossible parameter combinations          ####
    # #### Could possibly do this check in the loop at the function calls instead ####
    # if k_prime == 's2':
    #     if j > s2:
    #         return 0
    # else: # k_prime is s1
    #     if j > s1:
    #         return 0
    # ################################################################################

    W = (w_max + (s1 + s2)//2) - 1 # last position in seed

    if s1 <= s2:
        k_s = s1
        k_l = s2
    else:
        k_l = s1
        k_s = s2

    if (k_s != k_l) and (k_s > 0): # altstrobes
        if k_s == s1: # first strobe is short
            w_min = w_min - (k_l - k_s)//2
            w_max = w_max - (k_l - k_s)//2
        else: # first strobe is long
            w_min = w_min + (k_l - k_s)//2
            w_max = w_max + (k_l - k_s)//2


    # CASE 1 POSITION j IS ON THE FIRST STROBE
    if j <= s1:
        loc_i = 0 # first strobe
        if i <= (s1 - j): # i covered by first strobe
            p = 1
            loc_j = 0 # first strobe
        elif (s1 - j) < i < (w_min - j ): # within gap between first and second strobe
            p = 0
            loc_j = -1 # not cov
        elif (w_min - j ) <= i <= (W - j): # i within window of second strobe
            p = min(s2, (w_max - w_min + 1), i - (w_min - j) + 1, (W - j) - i + 1) / (w_max - w_min + 1)
            loc_j = 1 # second strobe
        else: # Outside valid range
            p = 0
            loc_j = -1 # not cov

    # CASE 2 POSITION j IS ON THE SECOND STROBE
    else:
        loc_i = 1 # second strobe
        if i <= s2 - (j - s1): # i covered by second strobe
            p = 1
            loc_j = 1 # second strobe
        else: # Outside second strobe
            p = 0
            loc_j = -1 # not cov

    # print('w_min=', w_min, 'W=', W, 'i=', i, 'j=', j, 'p=', p)
    return p, (loc_i, loc_j)



def compute_entropy_X_givenY(args, outfile):

    k = args.k
    w_min = args.w_min
    w_max = args.w_max
    W = w_max
    N = (w_max + k//2 - 2) # x range from 1 w-1 (w is max seed length in paper)


    ## GENERERAL THEORETICAL UPPER BOUND CALCULATION  NO LONGER WORKS UNDER OUR CONDITIONAL ENTROPY X|Y
    ##  because given a position j on the seed, we don't know the start position of the seed so we cannot
    ## provide a probability that remaining k-j+1 positions cover the reamingn positions in the window, as
    ## we don't know how many remaining positions in the window there are. This can be computed by also conditioning on
    ## Remaining positions. But we'll leave it to future work..

    # # Theoretical optimum (incorrect below under X|Y model)
    # X = [0]* N # maximum length of a fuzzy seed in this analysis
    # Xi_given_Yj = defaultdict(lambda: defaultdict(int))
    # for i in range(1, len(X)+1):
    #     for j in range(1,k+1):
    #         p_  = (k-j) / (w_max + k//2-1)  # k-j positions in the seed evenly distributed among  w positions
    #         Xi_given_Yj[j][i-1] += p_
    #     # X[j-1] = p_Xj # because 1-indexed
    # print(sum(X))
    # print('Optimal:', [round(x, 3) for x in X])
    # # print(Xi_given_Yj)
    # # approx from https://math.stackexchange.com/questions/2680787/entropy-of-a-binary-random-vector-with-exactly-k-non-zeros
    # # basically, entropy of X can be approximated by entropy of each individual X_i with a small correction factor. We have slightly different scenario though
    # sum_hx = 0
    # for s in Xi_given_Yj:
    #     X_c = Xi_given_Yj[s]
    #     for i in range(0, len(X_c)):
    #         if X_c[i] > 0 and (1- X_c[i]) > 0:
    #             sum_hx += X_c[i]* math.log(X_c[i],2) * (1/(k)) + (1- X_c[i]) * math.log((1- X_c[i]),2) * (1/(k))

    # H_X_theory_optimal  = -sum_hx
    # print("Optimal entropy: {0}".format( H_X_theory_optimal ) )

    # X2 = [(k) / (w_max + k//2-1)]*(w_max + k//2-1)
    # sum_hx = 0
    # for i in range(len(X2)):
    #     sum_hx += X2[i]* math.log(X2[i],2) + (1- X2[i]) * math.log((1- X2[i]),2)
    # print("Optimal entropy: {0}".format( -sum_hx ) )



    X = [0]* N # maximum length of a fuzzy seed in this analysis

    # k-mers

    Xi_given_Yj = defaultdict(lambda: defaultdict(int))
    for i in range(1, len(X)+1):
        for j in range(1,k+1):
            p_, (loc_i, loc_j) = P_X_i_given_Y_j(i, j, k, 0, w_min, w_max) # j is on k-mer
            Xi_given_Yj[j][i-1] += p_

    print(sum(X))
    print('K-mers:', [round(x, 3) for x in X])
    # print(Xi_given_Yj)
    # approx from https://math.stackexchange.com/questions/2680787/entropy-of-a-binary-random-vector-with-exactly-k-non-zeros
    # basically, entropy of X can be approximated by entropy of each individual X_i with a small correction factor. We have slightly different scenario though
    sum_hx = 0
    for s in Xi_given_Yj:
        X_c = Xi_given_Yj[s]
        for i in range(0, len(X_c)):
            if X_c[i] > 0 and (1- X_c[i]) > 0:
                sum_hx += X_c[i]* math.log(X_c[i],2) * (1/(k)) + (1- X_c[i]) * math.log((1- X_c[i]),2) * (1/(k))

    H_X_kmers  = -sum_hx

    print("k-mers entropy: {0}".format( H_X_kmers ) )
    outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'k-mers', k, '-', H_X_kmers, '-', 'altstrobes-analysis', 'Entropy'))


    # Altstrobes
    print("ALTSTROBES")
    for s in range(1, k//2+1):
        k_s = s
        k_l = k-s
        Xi_given_Yj = defaultdict(lambda: defaultdict(int))
        #print(s, w_min_short, w_max_short, w_min_long, w_max_long, (w_max_short - w_min_short), (w_max_long - w_min_long))
        for j in range(1,k+1):
            for i in range(1, len(X)+1):
                # p_Xj = 0
                p_, (loc_i, loc_j) = P_X_i_given_Y_j(i, j, k_s, k_l, w_min, w_max)
                if loc_i == 0 and loc_j == 0:
                    Xi_given_Yj[j][i-1] +=  0.5 * p_
                elif loc_i == 0 and loc_j == 1:
                    Xi_given_Yj[j][i-1] +=  0.5 * p_
                elif loc_i == 1 and loc_j == 1:
                    Xi_given_Yj[j][i-1] +=  0.5 * p_

                p_, (loc_i, loc_j) = P_X_i_given_Y_j(i, j, k_l, k_s, w_min, w_max)
                if loc_i == 0 and loc_j == 0:
                    Xi_given_Yj[j][i-1] +=  0.5 * p_
                elif loc_i == 0 and loc_j == 1:
                    Xi_given_Yj[j][i-1] +=  0.5 * p_
                elif loc_i == 1 and loc_j == 1:
                    Xi_given_Yj[j][i-1] +=  0.5 * p_

        sum_hx = 0
        for s in Xi_given_Yj:
            X_c = Xi_given_Yj[s]
            for i in X_c:
                if X_c[i] > 0 and (1- X_c[i]) > 0:
                    sum_hx += X_c[i]* math.log(X_c[i],2) * (1/(k)) + (1- X_c[i]) * math.log((1- X_c[i]),2) * (1/(k))

        H_X_altstrobes  = -sum_hx

        print("Altstrobes Entropy: {0}, ({1},{2})".format( H_X_altstrobes, k_s, k_l) )
        print()
        if k_s == k_l:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', k_s, '-', H_X_altstrobes, '-', 'altstrobes-analysis', 'Entropy'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'altstrobes', k_s, '-', H_X_altstrobes, '-', 'altstrobes-analysis', 'Entropy'))


    # Mixedstrobes
    for f in range(0, 11):
        p_strobe = f/10
        H_X = 0
        k_s = k//2
        k_l = k//2
        Xi_given_Yj = defaultdict(lambda: defaultdict(int))
        for i in range(1, len(X)+1):
            # p_Xj = 0
            for j in range(1,k+1):
                for t in ['kmer', 'strobemer']:
                    if t == 'kmer':
                        p_, (loc_i, loc_j) = P_X_i_given_Y_j(i, j, 0, k, w_min, w_max)
                        Xi_given_Yj[j][i-1] += ((1-p_strobe) * p_)
                    else:
                        p_, (loc_i, loc_j) = P_X_i_given_Y_j(i, j, k_s, k_l, w_min, w_max)
                        if loc_i == 0 and loc_j == 0:
                            Xi_given_Yj[j][i-1] += (p_strobe * p_)
                        elif loc_i == 0 and loc_j == 1:
                            Xi_given_Yj[j][i-1] += (p_strobe * p_)
                        elif loc_i == 1 and loc_j == 1:
                            Xi_given_Yj[j][i-1] += (p_strobe * p_)

        sum_hx = 0
        for s in Xi_given_Yj:
            X_c = Xi_given_Yj[s]
            for i in X_c:
                if X_c[i] > 0 and (1- X_c[i]) > 0:
                    sum_hx += X_c[i]* math.log(X_c[i],2) * (1/(k)) + (1- X_c[i]) * math.log((1- X_c[i]),2) * (1/(k))

        H_X_mixedstrobes  = -sum_hx

        print("Mixedstrobes Entropy: {0}, (frac:{1})".format( H_X_mixedstrobes, p_strobe) )
        print()

        if p_strobe == 0:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'k-mers', p_strobe, '-', H_X_mixedstrobes, '-', 'mixedstrobes-analysis', 'Entropy'))
        elif p_strobe == 1.0:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', p_strobe, '-', H_X_mixedstrobes, '-', 'mixedstrobes-analysis', 'Entropy'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'mixedstrobes', p_strobe, '-', H_X_mixedstrobes, '-', 'mixedstrobes-analysis', 'Entropy'))


    # Multistrobes
    print("MULTISTORBES")
    for k_boundary in range(1, k//2+1):
        Xi_given_Yj = defaultdict(lambda: defaultdict(int))
        p_k = 1/((k//2+1-k_boundary)) # p sampling a given length (k_s,k_l)
        print("k_boundary:", k_boundary, "p_k:", p_k)
        for i in range(1, len(X)+1):
            for j in range(1,k+1):
                for s in range(k_boundary, k//2+1): # possible strobe sizes
                    k_s = s
                    k_l = k-s

                    p_, (loc_i, loc_j) = P_X_i_given_Y_j(i, j, k_s, k_l, w_min, w_max)
                    if loc_i == 0 and loc_j == 0:
                        Xi_given_Yj[j][i-1] += p_k * 0.5 * p_
                    elif loc_i == 0 and loc_j == 1:
                        Xi_given_Yj[j][i-1] += p_k * 0.5 * p_
                    elif loc_i == 1 and loc_j == 1:
                        Xi_given_Yj[j][i-1] += p_k * 0.5 * p_

                    p_, (loc_i, loc_j) = P_X_i_given_Y_j(i, j, k_l, k_s, w_min, w_max)
                    if loc_i == 0 and loc_j == 0:
                        Xi_given_Yj[j][i-1] += p_k * 0.5 * p_
                    elif loc_i == 0 and loc_j == 1:
                        Xi_given_Yj[j][i-1] += p_k * 0.5 * p_
                    elif loc_i == 1 and loc_j == 1:
                        Xi_given_Yj[j][i-1] += p_k * 0.5 * p_

        sum_hx = 0
        for s in Xi_given_Yj:
            X_c = Xi_given_Yj[s]
            for i in X_c:
                if X_c[i] > 0 and (1- X_c[i]) > 0:
                    sum_hx += X_c[i]* math.log(X_c[i],2) * (1/(k)) + (1- X_c[i]) * math.log((1- X_c[i]),2) * (1/(k))

        H_X_multistrobes  = -sum_hx

        print("Multistrobes Entropy: {0}, (k_boundary: {1})\n".format( H_X_multistrobes, k_boundary) )

        if k_boundary == k//2:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', k_boundary, '-', H_X_multistrobes, '-', 'multistrobes-analysis', 'Entropy'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'continous altstrobes', k_boundary, '-', H_X_multistrobes, '-', 'multistrobes-analysis', 'Entropy'))



def main(args):
    outfile = open(args.outcsv, 'w')
    outfile.write("w_min,w_max,type,x,m,y,mean_error_y_est,analysis,ymetric\n") # ymetric either p_match, E_overlap_theory, E_overlap_real_hash, Entropy
    print('Using N_SIM:', args.N_SIM)
    compute_entropy_X_givenY(args, outfile)
    compute_p_match(args, args.N_SIM, outfile)
    outfile.close()

    outfile = open(args.outcsv, 'r')
    plot_p_match(args, outfile.name, args.outplot_prefix)
    plot_entropy(args, outfile.name, args.outplot_prefix)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    parser.add_argument('k', type=int, default=15, help='k-mer size')
    parser.add_argument('w_min', type=int, default=25, help='Window size')
    parser.add_argument('w_max', type=int, default=50, help='Window size')
    # parser.add_argument('f', type=str, default="2/3", help='relative size between strobes within strobemer')
    parser.add_argument('outcsv', type=str,  default=None, help='output CSV file for plotting.')
    parser.add_argument('outplot_prefix', type=str,  default=None, help='output pdf file prefix.')
    parser.add_argument('--max_mut_freq', type=float, default=0.3, help='Maximum mutation frequency for simulations')
    parser.add_argument('--N_SIM', type=int, default=1000, help='Number of simulations for p_matches experiment')

    # parser.set_defaults(which='main')
    args = parser.parse_args()



    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit()

    main(args)
