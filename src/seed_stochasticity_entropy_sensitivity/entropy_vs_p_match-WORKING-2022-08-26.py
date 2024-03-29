import os,sys
import argparse

import random
import math
from collections import defaultdict

import indexing_Maier_altstrobes as indexing
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot


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

    g = sns.relplot(data=indata_entropy, x="x", y="y", hue="type", col="analysis", s=80, facet_kws=dict(sharex=False))
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

    g = sns.relplot(data=indata_p_match, x="x", y="y", hue="type", style="m", col="analysis", s=80, facet_kws=dict(sharex=False))
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
    
    N_m = [i for i in range(1,61)] 
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

    ###################################


    for m in N_m:
        print('m=',m)
        seq_pairs = []
        for n in range(N_sim):
            seq1 = "".join([random.choice("ACGT") for i in range(L)])
            muts = set(random.sample(range(2*W),m)) # only sample mutations within the 2W region to control mutation rate in the region according to theory 
            seq2 = "".join([seq1[i] if i not in muts else random.choice(['', reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
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
                altstrobes_pos1 = indexing.altstrobes(seq1, k_s, k_l, w_min, w_max, 1)
                # print(altstrobes_pos1)
                altstrobes_pos2 = indexing.altstrobes(seq2, k_s, k_l, w_min, w_max, 1)
                
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
            print("mixedstrobes-{0}: p_match: {1}, m:{2}".format(p_strobe,2, p_match, m) )

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


# def P_X_j_given_X_i_and_c_i(j, c_i, k_s, k_l, w_min, w_max, k_prime = 'k_l'):
#     """
#         k_prime indicates which strobe is currnetly covering position i
#         j is given as offset from i which is assumed to be at position 1. 

#         Assumed j > i here.
#         c_i is in [1,k]
#         Also assumes w_min > k_l
#     """

#     ################################################################################
#     #### Below if-else takes care of impossible parameter combinations          ####
#     #### Could possibly do this check in the loop at the function calls instead #### 
#     if k_prime == 'k_l':
#         if c_i > k_l:
#             return 0
#     else: # k_prime is k_s
#         if c_i > k_s:
#             return 0
#     ################################################################################

#     W = (w_max + (k_s + k_l)//2) # last position in seed

#     if k_prime == 'k_l':
#         if k_l >= j: # j covered by k_l, i is 1 indexed
#             p = 1
#         elif w_min <= j: # j covered by k_s i is 1 indexed
#             p = min(k_s, j - w_min + 1, W - j) / (w_max - w_min) 
#         else: # within gap between first and second strobe
#             p = 0
#         # print('here',k_l,c_i, j, p)
#     else:
#         if k_s >= j: # j covered by k_l, i is 1 indexed
#             p = 1
#         elif w_min <= j: # j covered by k_s i is 1 indexed
#             p = min(k_l, j - w_min + 1, W - j) / (w_max - w_min) 
#         else: # within gap between first and second strobe
#             p = 0

#     return p


def P_X_j_given_X_i_and_c_i_NEW(j, c_i, s1, s2, w_min, w_max):
    """
        s1: strobe 1 length.
        s2: strobe 2 length.

        k_prime indicates which strobe is currently covering position i
        j is given as offset from i which is assumed to be at position 1. 

        Assumed j > i here.
        c_i is in [1,k]
        Also assumes w_min > s2
    """

    # ################################################################################
    # #### Below if-else takes care of impossible parameter combinations          ####
    # #### Could possibly do this check in the loop at the function calls instead #### 
    # if k_prime == 's2':
    #     if c_i > s2:
    #         return 0
    # else: # k_prime is s1
    #     if c_i > s1:
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


    # CASE 1 POSITION C_I IS ON THE FIRST STROBE
    if c_i <= s1:
        if j <= (s1 - c_i): # j covered by first strobe
            p = 1
        elif (s1 - c_i) < j < (w_min - (c_i - 1)): # within gap between first and second strobe
            p = 0
        elif (w_min - (c_i - 1)) <= j <= (W - c_i): # j covered by second strobe 
            # print(s2, j, c_i, 'min({0}, {1}, {2})'.format(s2, j - (w_min - c_i) + 1, (W - c_i) - j + 1))
            p = min(s2, (w_max - w_min), j - (w_min - c_i) + 1, (W - c_i) - j + 1) / (w_max - w_min) 
        else: # Outside valid range
            p = 0

    # CASE 2 POSITION C_I IS ON THE SECOND STROBE
    else:
        if j <= s2 - (c_i - s1): # j covered by second strobe
            p = 1
        else: # Outside second strobe
            p = 0   

    # print('w_min=', w_min, 'W=', W, 'j=', j, 'c_i=', c_i, 'p=', p)
    return p



def compute_entropy_Xj(args, outfile):
    """
        # Assumes non-overlapping strobes for all strobemers, i.e., k unique positions sampled

        X_i = Event that position i are covered by a seed (Binary variable)
        P( X_j | X_i ) = Probability that position j are covered by the seed if i is covered
        Let c_i be the coordinate within the seed that covers i.

        P( X_j ) &= \sum_{c = 0}^{k-1} \sum_{m = 1}^M P( X_j | Y_m, c) P(Y_m, c) \\
                       &=  \sum_{c = 0}^{k-1} \sum_{m = 1}^M P( X_j | Y_m, c)  P( Y_m | c) P(c) \\
                      &=   \frac{1}{k} \sum_{c_i = 0}^{k-1} \sum_{m = 1}^M P( X_j | Y_m, c)  P( Y_m | c).


        Above holds for k-mers, altstrobes and randstrobes. For mixedstrobes we need to condition on seed type. We have
        P( X_j | X_i ) = \sum_{c_i =1}^{30} sum_{k_i \in k_s,k_l} sum_{t \in k-mer, strobemer} P( X_j | X_i, c_i, k_i, t) * P(c_i, k_i, t) = 
                         \sum_{c_i =1}^{30} sum_{k_i \in k_s,k_l} sum_{t \in k-mer, strobemer} P( X_j | X_i, c_i, k_i, t) * P(c_i | k_i, t) * P(k_i, t)               =
                         \sum_{c_i =1}^{30} sum_{k_i \in k_s,k_l} sum_{t \in k-mer, strobemer} P( X_j | X_i, c_i, k_i, t) * P(c_i | k_i, t) * P(k_i | t)       * p(t) =
                         \sum_{c_i =1}^{30} sum_{k_i \in k_s,k_l} sum_{t \in k-mer, strobemer} P( X_j | X_i, c_i, k_i, t) * (1/k_i)         * (k_i/(k_s + k_l) * p(t) =
                         \sum_{c_i =1}^{30} sum_{k_i \in k_s,k_l} sum_{t \in k-mer, strobemer} P( X_j | X_i, c_i, k_i, t) * (1/(k_s + k_l) * p(t)

        From above, we can compute the entropy of X_j | X_i exactly.
        However, for X = \sum_{j=0}^{|S|} = \sum_{j=i-W-k}^{i+W+k} = X_j there is no feasable calculation. And furthermore, the X_j are not independent and exactly k  X_j can be 1 at a time.
        We rely on the following approximation: https://math.stackexchange.com/questions/2680787/entropy-of-a-binary-random-vector-with-exactly-k-non-zeros

    """


    k = args.k
    w_min = args.w_min
    w_max = args.w_max
    W = w_max


    # Theoretical optimum
    X = [ k/(w_max + k//2) ] * (w_max + k//2) # maximum length of a fuzzy seed in this analysis
    print('Optimum:', [round(x, 3) for x in X])
    X2 = [ 0.5 ] * (w_max + k//2) # maximum length of a fuzzy seed in this analysis


    # approx from https://math.stackexchange.com/questions/2680787/entropy-of-a-binary-random-vector-with-exactly-k-non-zeros
    # basically, entropy of X can be approximated by entropy of each individual X_i with a small correction factor. We have slightly different scenario though
    sum_hx = 0
    sum_hx2 = 0
    for i in range(0, len(X)):
        if X[i] > 0 and (1- X[i]) > 0:
            sum_hx += X[i]* math.log(X[i],2) + (1- X[i]) * math.log((1- X[i]),2)
            sum_hx2 += X2[i]* math.log(X2[i],2) + (1- X2[i]) * math.log((1- X2[i]),2)
    
    H_X_optimum  = -sum_hx
    H_X_optimum2  = -sum_hx2
    
    print("Optimal entropy: {0}".format( H_X_optimum ) )
    print("Optimal entropy general (no restrct on k): {0}".format( H_X_optimum2 ) )
    for k_s in range(1, k+1):
        outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'Optimal', k_s, '-', H_X_optimum, '-', 'altstrobes-analysis', 'Entropy'))

    X = [0]* (w_max + k//2) # maximum length of a fuzzy seed in this analysis
    
    # k-mers
    j = w_max + k
    for j in range(1, len(X)):
        p_Xj = 0
        for c_i in range(1,k+1):
            # print(j, c_i, P_X_j_given_X_i_and_c_i_NEW(j, c_i, 0, k, w_min, w_max, k_prime = 'k_l'))
            p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, k, 0, w_min, w_max) * (1/(k)) # c_i is on k-mer
            #p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, 0, k, w_min, w_max, k_prime = 'k_s') * (1/(k)) # c_i is on short strobe
        X[j-1] = p_Xj # because 1-indexed

    print('K-mers:', [round(x, 3) for x in X])

    # approx from https://math.stackexchange.com/questions/2680787/entropy-of-a-binary-random-vector-with-exactly-k-non-zeros
    # basically, entropy of X can be approximated by entropy of each individual X_i with a small correction factor. We have slightly different scenario though
    sum_hx = 0
    for i in range(0, len(X)):
        if X[i] > 0 and (1- X[i]) > 0:
            sum_hx += X[i]* math.log(X[i],2) + (1- X[i]) * math.log((1- X[i]),2)
    
    H_X_kmers  = -sum_hx
    
    print("k-mers entropy: {0}".format( H_X_kmers ) )
    outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'k-mers', k, '-', H_X_kmers, '-', 'altstrobes-analysis', 'Entropy'))


    # Altstrobes
    print("ALTSTROBES")
    X = [0]* (w_max + k//2) # maximum length of a fuzzy seed in this analysis
    for s in range(1, k//2+1):
        k_s = s
        k_l = k-s
        sum_hx = 0
        #print(s, w_min_short, w_max_short, w_min_long, w_max_long, (w_max_short - w_min_short), (w_max_long - w_min_long))
        for j in range(1, len(X)):
            p_Xj = 0
            for c_i in range(1,k+1):
                # print(j, c_i, P_X_j_given_X_i_and_c_i_NEW(j, c_i, 0, k, w_min, w_max, k_prime = 'k_l'))
                # p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_s, k_l, w_min_short, w_max_short, k_prime = 'k_l') * (1/(k)) # c_i is on long strobe
                # p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_s, k_l, w_min_long, w_max_long, k_prime = 'k_s') * (1/(k)) # c_i is on short strobe
                p_Xj += 0.5*P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_s, k_l, w_min, w_max) * (1/k)
                p_Xj += 0.5*P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_l, k_s, w_min, w_max) * (1/k)
                # p_Xj += 0.5 * P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_s, k_l, w_min_short, w_max_short) * (1/(k)) # first strobe is short
                # p_Xj += 0.5 * P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_l, k_s, w_min_long, w_max_long) * (1/(k)) # first strobe is long
            X[j-1] = p_Xj # because 1-indexed

        print('Altstrobes:' , [round(x, 3) for x in X])
        #print('0-positions:' , [i for i, x in enumerate(X) if x == 0])
        #sys.exit()
        
        sum_hx = 0
        for i in range(0, len(X)):
            if X[i] > 0 and (1- X[i]) > 0:
                sum_hx += X[i]* math.log(X[i],2) + (1- X[i]) * math.log((1- X[i]),2)

        H_X_altstrobes  = -sum_hx
        
        print("Altstrobes Entropy: {0}, ({1},{2})".format( H_X_altstrobes, k_s, k_l) )
        print()
        if k_s == k_l:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', k_s, '-', H_X_altstrobes, '-', 'altstrobes-analysis', 'Entropy'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'altstrobes', k_s, '-', H_X_altstrobes, '-', 'altstrobes-analysis', 'Entropy'))


    # Mixedstrobes
    X = [0]* (w_max + k//2) # maximum length of a fuzzy seed in this analysis
    for f in range(0, 11):
        p_strobe = f/10
        H_X = 0
        k_s = k//2
        k_l = k//2
        for j in range(1, len(X)):
            p_Xj = 0
            for c_i in range(1,k+1):
                for t in ['kmer', 'strobemer']:
                    if t == 'kmer':
                        # p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, 0, k, w_min, w_max, k_prime = 'k_l') * (1/k) * (1-p_strobe) # c_i is on long strobe
                        # p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, 0, k, w_min, w_max, k_prime = 'k_s') * (1/k) * (1-p_strobe) # c_i is on short strobe
                        p_Xj += (1-p_strobe) * P_X_j_given_X_i_and_c_i_NEW(j, c_i, 0, k, w_min, w_max) * (1/k) # c_i is on k-mer
                    else:
                        # p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_s, k_l, w_min, w_max, k_prime = 'k_l') * (1/(k)) * p_strobe # c_i is on long strobe
                        # p_Xj += P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_s, k_l, w_min, w_max, k_prime = 'k_s') * (1/(k)) * p_strobe # c_i is on short strobe
                        p_Xj +=  p_strobe * P_X_j_given_X_i_and_c_i_NEW(j, c_i, k_s, k_l, w_min, w_max) * (1/(k))

            X[j-1] = p_Xj # because 1-indexed

        print('Mixedstrobes:' , [round(x, 3) for x in X])
        
        sum_hx = 0
        for i in range(0, len(X)):
            if X[i] > 0 and (1- X[i]) > 0:
                sum_hx += X[i]* math.log(X[i],2) + (1- X[i]) * math.log((1- X[i]),2)

        H_X_mixedstrobes  = -sum_hx
        
        print("Mixedstrobes Entropy: {0}, (frac:{1})".format( H_X_mixedstrobes, p_strobe) )
        print()

        if p_strobe == 0:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'k-mers', p_strobe, '-', H_X_mixedstrobes, '-', 'mixedstrobes-analysis', 'Entropy'))
        elif p_strobe == 1.0:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'randstrobes', p_strobe, '-', H_X_mixedstrobes, '-', 'mixedstrobes-analysis', 'Entropy'))
        else:
            outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(w_min, w_max, 'mixedstrobes', p_strobe, '-', H_X_mixedstrobes, '-', 'mixedstrobes-analysis', 'Entropy'))



def main(args):
    outfile = open(args.outcsv, 'w')
    outfile.write("w_min,w_max,type,x,m,y,mean_error_y_est,analysis,ymetric\n") # ymetric either p_match, E_overlap_theory, E_overlap_real_hash, Entropy
    print('Using N_SIM:', args.N_SIM)
    compute_entropy_Xj(args, outfile)
    compute_p_match(args, args.N_SIM, outfile)
    outfile.close()
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
    parser.add_argument('--N_SIM', type=int, default=1000, help='Number of simulations for p_matches experiment')

    # parser.set_defaults(which='main')
    args = parser.parse_args()



    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit()

    main(args)