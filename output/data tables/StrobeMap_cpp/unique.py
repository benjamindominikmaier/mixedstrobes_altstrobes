method = "mixedminstrobe3_0.9"
f = open(method+".tsv", "r")

matches = dict()
for line in f:
    if line.strip().startswith(">"):
        id = line
        matches[id] = {}
    else:
        r_id, r_pos, q_pos, r_len = line.strip().split(" ")

        r_pos = int(r_pos)
        q_pos = int(q_pos)
        r_len = int(r_len)

        if not (r_pos+r_len, q_pos+r_len) in matches[id]:
            matches[id][(r_pos+r_len, q_pos+r_len)] = (r_id, r_pos, q_pos, r_len)
        elif r_len > matches[id][(r_pos+r_len, q_pos+r_len)][3]:
            matches[id][(r_pos+r_len, q_pos+r_len)] = (r_id, r_pos, q_pos, r_len)
        else:
            continue

outfile = open(method+"_fixed.tsv", 'w')
for id in matches.keys():
    outfile.write("{0}".format(id))
    for elem in list(matches[id].values()):
        ref_acc, ref_p, q_pos, k = elem
        outfile.write("  {0} {1} {2} {3}\n".format(ref_acc, ref_p, q_pos, k))
