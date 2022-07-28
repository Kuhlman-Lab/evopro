import json
import re
import random

class DesignSeq:

    def __init__(self, jsonfile=None, seq=None, mutable=None, sym=[]):
        aalists=[]
        resids=[]

        if jsonfile is not None:
            s, m = self.load_DesignSeq_json(jsonfile)
            self.seq = s
            self.mutable = m
            with open(jsonfile, "r") as inpf:
                self.jsondata = json.load(inpf)
        else:
            self.seq = seq
            self.mutable = mutable
            self.jsondata = self.create_jsondata(seq, mutable, sym)

        if self.mutable is not None:
            for chain in self.seq:
                for aa, i in zip(self.seq[chain], range(len(self.seq[chain]))):
                    resid = chain+str(i+1)
                    resids.append(resid)

                    if resid in self.mutable:
                        aalists.append(self.mutable[resid])
                    else:
                        aalists.append(None)
        else:
            print("Warning: mutable residues not specified. All residues will be designed")
            self.mutable={}
            for chain in self.seq:
                for aa, i in zip(self.seq[chain], range(len(self.seq[chain]))):
                    resid = chain+str(i+1)
                    resids.append(resid)
                    aalists.append([1 for x in range(20)])
                    self.mutable[resid] = [1 for x in range(20)]

        self.aalists = aalists
        self.resids = resids

    def __eq__(self, other):
        return self.seq == other.seq

    def __str__(self):
        return (str(self.seq) + "\n" + str(self.mutable) + "\n" + str(self.jsondata))

    def __hash__(self):
        return hash(str(self.resids) + "\t" + str(self.aalists))

    def load_DesignSeq_json(self, jsoninputsfile):
        all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        with open(jsoninputsfile, "r") as inpf:
            data = json.load(inpf)

        seq = data["sequence"]
        mut_res = {}
        #ADD FUNCTIONALITY FOR + RESIDUES FOR EACH OPTION?
        for res in data["designable"]:
            resid = res["chain"]+str(res["resid"])
            weighted_aas = []
            w = []
            if "hydphob" in res["MutTo"]:
                w = [1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1]
            elif "hydphil" in res["MutTo"]:
                w = [0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            elif "alpha" in res["MutTo"]:
                w = [1, 0, 0.31, 0.6, 0.46, 0, 0.39, 0.59, 0.74, 0.79, 0.74, 0.35, 0, 0.61, 0.79, 0.5, 0.34, 0.39, 0.51, 0.47]
            elif "all" in res["MutTo"]:
                w = [1 for x in all_aas]
            removed = []
            if "-" in res["MutTo"]:
                removed = list(res["MutTo"].split("-")[1])
            for aa, i in zip(all_aas, range(len(all_aas))):
                if aa in removed:
                    weighted_aas.append(0)
                else:
                    weighted_aas.append(w[i])

            if "custom" in res["MutTo"]:
                try:
                    weighted_aas = [float(x) for x in str(res["aalist"]).strip("[]").split(",")]
                except:
                    print("Custom weights not specified. Defaulting to all")
                    weighted_aas = [1 for x in range(20)]

            if not weighted_aas:
                weighted_aas = [1 for x in range(20)]
                print("User did not specify valid residue substitutions for residue " + res["chain"]+str(res["resid"]) + ". Defaulting to all")
             
            mut_res[resid] = weighted_aas

        return seq, mut_res

    def create_jsondata(self, seq, mutable, sym):
        json_dict = {"sequence":seq}
        des = []
        for resid in mutable:
            chain = re.split('(\d+)', resid)[0]
            aa = int(re.split('(\d+)', resid)[1])
            aa_id = self.get_aa_identity(resid)
            mutto = "custom"
            aalist = self.mutable[resid]
            des_dict = {"chain":chain, "resid":aa, "WTAA":aa_id, "MutTo":mutto, "aalist": aalist}
            des.append(des_dict)
        json_dict["designable"] = des
        json_dict["symmetric"] = sym
        return json.dumps(json_dict)

    def mutate(self, mutpercent = 0.125, symmetry=None):
        all_aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        #calculating number of mutants from mutpercent and maybe adding or subtracting one mutant for stochasticity
        #num_mut = round(len(str(self))*mutpercent) + random.choice([-1, 0, 1])
        num_mut=round(mutpercent*len(self.mutable.keys()))+random.choice([-1, 0, 1])
        if num_mut<1:
            num_mut=1

        oldseq = {chain: list(self.seq[chain]) for chain in self.seq}
        newseq = {chain: list(self.seq[chain]) for chain in self.seq}

        if symmetry is None:
            symmetry = []

        mut_ids = random.sample(list(self.mutable.keys()), num_mut)
        for mut_id in mut_ids:
            weights = self.mutable[mut_id]
            chain = re.split('(\d+)', mut_id)[0]
            aa = int(re.split('(\d+)', mut_id)[1])

            new_aa = random.choices(all_aas, weights)[0]
            newseq[chain][aa-1] = new_aa
            if chain in symmetry:
                for other_chain in symmetry:
                    newseq[other_chain] = newseq[chain]

        newseq_joined = {chain: "".join(newseq[chain]) for chain in newseq}
        newseqobj = DesignSeq(seq=newseq_joined, mutable=self.mutable)
        return newseqobj

    def get_aa_identity(self, index):
        chain = re.split('(\d+)', index)[0]
        aa = int(re.split('(\d+)', index)[1])
        return self.seq[chain][aa-1]

    def crossover(self, otherDS, ncross = 1):
        mut_seq = sorted(self.mutable.keys())
        seq1 = [self.get_aa_identity(mut_id) for mut_id in mut_seq]
        other_mut_seq = sorted(otherDS.mutable.keys())
        seq2 = [otherDS.get_aa_identity(mut_id) for mut_id in other_mut_seq]
        points = random.sample(range(len(seq1)), ncross)
        s=0
        newseq = []
        newseqforobj = {chain: list(self.seq[chain]) for chain in self.seq}
        for i, aa1, aa2 in zip(range(len(seq1)), seq1, seq2):
            options = [aa1, aa2]
            if i in points:
                if s==0:
                    s=1
                else:
                    s=0
            newseq.append(options[s])

        #print(seq1, seq2, newseq)

        for mut_id, new_aa in zip(mut_seq, newseq):
            chain = re.split('(\d+)', mut_id)[0]
            aa = int(re.split('(\d+)', mut_id)[1])
            newseqforobj[chain][aa-1] = new_aa

        newseq_joined = {chain: "".join(newseqforobj[chain]) for chain in newseqforobj}
        newseqobj = DesignSeq(seq=newseq_joined, mutable=self.mutable)
        return newseqobj

if __name__ == "__main__":
    dsobj = DesignSeq(jsonfile="/pine/scr/a/m/amritan/kuhlmanlab/fdd/fdd/tests/pd1_threehelix_test/resfile.json")
    dsobj2 = DesignSeq(jsonfile="test2.json")
    newdsobj = dsobj.mutate()
    newdsobj2 = dsobj.crossover(dsobj2)
    print(dsobj.seq, newdsobj.seq, newdsobj2.seq)
