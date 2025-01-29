oldfile = "inputs/sequences.csv"

newseqs = []
with open(oldfile, "r") as f:
    for lin in f:
        seq = lin.strip("\n,").split(",")
        newseq = []
        # newseq = newseq + [s + "QSFYDLSLQTSEIKEEEKELRRTTYKLQVKN" for s in seq[:3]]
        # newseq = newseq + seq[3:]
        newseq = newseq  + [seq[0]] + [seq[0]] + [seq[0]]
        #newseq = newseq + ["PEPKSRFAMLDDVKILANGLLQLGHGLKDFVHKTKGQINDIFQKLNIFDQSFYDLSLQTSEIKEEEKELRRTTYKLQVKN" for k in range(3)]
        newseq = newseq + ["GPVQSKSPRFASWDEMNVLAHGLLQLGQGLREHAERTRSQLSALERRLSACGSACQGTE" for k in range(3)]
        newseqstr = ",".join(newseq)
        print(seq,newseqstr)
        newseqs.append(newseqstr)

with open("newfile.csv", "w") as opf:
    for n in newseqs:
        opf.write("," + n + "\n")
