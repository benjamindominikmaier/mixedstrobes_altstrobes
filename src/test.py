from genome_mapping_metrics import *

read_coverage_solution = {}
total_disjoint_matches = 0
tot_genome_length = 0

for (query_acc, nams) in get_NAM_records("tmp2.tsv", {"NZ_CP027472.1": 1}):
    q_acc = query_acc.split()[0]
    tot_genome_length += 500
    # print(t_chr_id, t_start, t_end)
    for ref_id in nams:
        print("Ref ID", ref_id)
        chrom_nams = nams[ref_id]
        total_disjoint_matches += len(chrom_nams)
        solutions, opt_cov = n_logn_read_coverage(chrom_nams)
        #print("SOLUTION", solutions)
        # print("OPT_COV", opt_cov)
        # pick best from forward and reverse strand
        if q_acc in read_coverage_solution:
            c_prev, _ = read_coverage_solution[q_acc]
            if c_prev < opt_cov:
                # print(ref_id, opt_cov)
                read_coverage_solution[q_acc] = (opt_cov, solutions[0])
        else:
            # print(ref_id, opt_cov)
            read_coverage_solution[q_acc] = (opt_cov, solutions[0])


# print(read_coverage_solution)
tot_genome_length = tot_genome_length/2 # remove double counting of reverse complements
collinear_chain_nam_sizes = []
total_bp_covered = 0

sc_positions = []
mc_positions = []

# collinear_outfile = open("collinear_matches_out.tsv", "w")
gap_pos = 0
gap_length = 0
gaps = []
for q_acc in read_coverage_solution:
    # collinear_outfile.write("> {0}\n".format(q_acc))
    opt_cov, solution = read_coverage_solution[q_acc]
    total_bp_covered += opt_cov
    for n in solution:
        collinear_chain_nam_sizes.append(n.y - n.x)
        sc_positions.append(n.c)
        sc_positions.append(n.c + n.val - 10)

        if n.c > gap_pos:
            print(n.c, n.d)
            gaps.append(n.c-gap_pos-1)
            gap_length += n.c-gap_pos-1
            print("GAP", n.c-gap_pos-1)
        gap_pos = n.d
        # collinear_outfile.write("  {0} {1} {2} {3}\n".format(n.chr_id, n.x, n.c, n.val))

gap_length += 500 - gap_pos
print(gap_length)
print(total_bp_covered)
coll_esize = e_size(collinear_chain_nam_sizes, tot_genome_length)
# collinear_outfile.close()
