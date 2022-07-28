

def write_starting_seqs(filename, seqslist):
    with open(filename, "w") as opf:
        for seq in seqslist:
            opf.write(seq.strip() + "\n")

def read_seqs_finalscores(finalscoresfile):
    seqslist = []
    with open(finalscoresfile, "r") as f:
        for lin in f:
            l = lin.strip().split("\t")
            print(l)
            seqslist.append(l[1])
    return seqslist

if __name__ == "__main__":
    seqs = read_seqs_finalscores("/pine/scr/a/m/amritan/kuhlmanlab/fdd/fdd/tests/pd1_threehelix_run1/final_seqs_scores.txt")
    write_starting_seqs("/pine/scr/a/m/amritan/kuhlmanlab/fdd/fdd/tests/pd1_threehelix_run3/starting_seqs.txt", seqs)
