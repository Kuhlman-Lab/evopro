import random

one_letter = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

three_letter = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", 
                "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

def three_to_one(code):
    for one, three in zip(one_letter, three_letter):
        if code == three:
            return one
    
def one_to_three(code):
    for one, three in zip(one_letter, three_letter):
        if code == one:
            return three

def generate_randoms(num_seqs, length):
    seqs = []
    for i in range(num_seqs):
        seqs.append("".join(random.choices(one_letter, k=length)))

    return seqs

if __name__ == "__main__":
    #print(three_to_one("ASN"))
    #print(one_to_three("S"))
    seqs = generate_randoms(19, 59)
    for seq in seqs:
        print(seq)
