oldfile = "inputs/sequences.csv"

newseqs = []
with open(oldfile, "r") as f:
    for lin in f:
        seq = lin.strip("\n,").split(",")
        newseq = []
        newseq = newseq + [s + "QSFYDLSLQTSEIKEEEKELRRTTYKLQVKN" for s in seq[:3]]
        newseq = newseq + seq[3:]
        newseqstr = ",".join(newseq)
        print(seq,newseqstr)
        newseqs.append(newseqstr)

with open("newfile.csv", "w") as opf:
    for n in newseqs:
        opf.write("," + n + "\n")
