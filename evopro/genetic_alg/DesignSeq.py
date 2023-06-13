import json
import re
import random
import copy

class DesignSeq:

    def __init__(self, jsonfile=None, sequence=None, mutable=None, symmetric=[], jdata=None, seq=None):

        if jsonfile is not None:
            jseq, mut, sym, jdata = self._load_from_json(jsonfile)
            self.sequence = jseq
            self.mutable = mut
            self.symmetric = sym
            self.jsondata = jdata
            self._create_jsondata()

        elif seq is not None:
            self.sequence = sequence
            self.mutable = mutable
            self.mutable = self._update_mutable_from_seq(seq)
            self._update_sequence()
            self.symmetric = symmetric
            self._create_jsondata()
        else:
            self.sequence = sequence
            self.mutable = mutable
            self.symmetric = symmetric
            self._update_sequence()
            self._create_jsondata()

    def __eq__(self, other):
        return self.get_sequence_string() == other.get_sequence_string()

    def __str__(self):
        return self.get_sequence_string()

    def __hash__(self):
        return hash(self.get_sequence_string())

    def _load_mutable(self, designable):
        """takes user input of the designable residue json list and beefs up internal representation"""

        all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        d = copy.deepcopy(designable)
        mut = {}

        for res in d:
            res["resid"] = [res["chain"], res["resid"], 0, res["WTAA"]]
            weighted_aas = []
            w = []
            if "hydphob" in res["MutTo"]:
                w = [1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1]
            elif "hydphil" in res["MutTo"]:
                w = [0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            elif "alpha" in res["MutTo"]:
                w = [1, 0, 0.31, 0.6, 0.46, 0, 0.39, 0.59, 0.74, 0.79, 0.74, 0.35, 0, 0.61, 0.79, 0.5, 0.34, 0.39, 0.51, 0.47]
            elif "all" in res["MutTo"]:
                w = [1 for x in all_aas]
            elif set(res['MutTo']).issubset(set('ACDEFGHIKLMNPQRSTVWY')):
                w = [0 for x in all_aas]
                for aa in set(res['MutTo']):
                    w[all_aas.index(aa)] = 1

            removed = []
            if "-" in res["MutTo"]:
                removed = list(res["MutTo"].split("-")[1])
            added = []
            if "+" in res["MutTo"]:
                added = list(res["MutTo"].split("+")[1])
            for aa, i in zip(all_aas, range(len(all_aas))):
                if aa in added:
                    weighted_aas.append(1)
                elif aa in removed:
                    weighted_aas.append(0)
                else:
                    weighted_aas.append(w[i])

            """
            #code for handling custom weights...needs rewrite for new rep
            #maybe have a separate file for specifying custom weights?
            if "custom" in res["MutTo"]:
                try:
                    weighted_aas = [float(x) for x in str(res["aalist"]).strip("[]").split(",")]
                except:
                    print("Custom weights not specified. Defaulting to all")
                    weighted_aas = [1 for x in range(20)]
            """

            if not weighted_aas:
                weighted_aas = [1 for x in range(20)]
                print("User did not specify valid residue substitution options for residue " + res["chain"]+str(res["resid"]) + ". Defaulting to all")
            
            del res["WTAA"]
            del res["chain"]
            res["weights"] = weighted_aas
            key = res["resid"][0] + str(res["resid"][1])
            mut[key] = [res]

        return mut

    def _load_sequence(self, seq_dict):
        """takes user input of the sequence from json and beefs up internal representation"""
        seq = {}
        for chain in seq_dict:
            sequence = seq_dict[chain]
            for pos, res in zip(range(len(sequence)), sequence):
                resid = chain + str(pos+1)
                seq[resid] = res

        return seq

    def _load_symmetry(self, sym_list):
        """takes user input of the symmetry list from json and beefs up internal representation"""
        sym_dict = {}
        for slist in sym_list:
            for elem1 in slist:
                if elem1 not in sym_dict:
                    sym_dict[elem1] = []
                for elem2 in slist:
                    if elem1 != elem2 and elem2 not in sym_dict[elem1]:
                        sym_dict[elem1].append(elem2)
        return sym_dict

    def _load_from_json(self, jsoninputsfile):
        """load a DesignSeq object from a json file with specifications"""
        with open(jsoninputsfile, "r") as inpf:
            jdata = json.load(inpf)

        if "symmetric" not in jdata:
            jdata["symmetric"] = []

        seq = self._load_sequence(jdata["sequence"])
        mut = self._load_mutable(jdata["designable"])
        sym = self._load_symmetry(jdata["symmetric"])


        return seq, mut, sym, jdata

    def _create_jsondata(self):
        """creates json dictionary for dsobj that can be used as input to protein mpnn"""
        all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        new_des = []
        new_sym = []
        new_jsonseq = {}

        for resid in self.sequence:
            chain = re.split('(\d+)', resid)[0]
            if chain not in new_jsonseq:
                new_jsonseq[chain] = []

            new_jsonseq[chain].append(self.sequence[resid])

        newjsonseq_joined = {chain: "".join(new_jsonseq[chain]) for chain in new_jsonseq}
        json_dict = {"sequence":newjsonseq_joined}

        chains = []
        numbering = {}
        i = 1
        for resid in self.sequence:
            chain = re.split('(\d+)', resid)[0]
            #checking if we need to restart the numbering from 1 because new chain
            if chain not in chains:
                chains.append(chain)
                i = 1
                numbering[chain] = []

            resval = self.sequence[resid]
            if resid in self.mutable:
                for res in self.mutable[resid]:
                    if res["resid"][3] != '':
                        des_dict = {"chain":chain, "resid":i, "WTAA":res["resid"][3], "MutTo":res["MutTo"]}
                        new_des.append(des_dict)
                        #print("designable", des_dict)
                        numbering[chain].append(resid)
                        i+=1

                    else:
                        pass
            else:
                numbering[chain].append(resid) 
                i+=1


        json_dict["designable"] = new_des

        sym_sets = []
        chains = []
        inc = 0
        for resid in self.sequence:
            if resid in self.symmetric:
                chain1 = re.split('(\d+)', resid)[0]
                resnum1 = int(re.split('(\d+)', resid)[1])
                if len(self.sequence[resid]) == 1:
                    present = False
                    for s in sym_sets:
                        if chain1+str(resnum1+inc) in s or (resnum1+inc)<1:
                            present = True
                            break
                    if not present:
                        new_set = {chain1+str(resnum1+inc)}
                        for elem in self.symmetric[resid]:
                            chain2 = re.split('(\d+)', elem)[0]
                            resnum2 = int(re.split('(\d+)', elem)[1])
                            new_set.add(chain2+str(resnum2+inc))
                        sym_sets.append(new_set)
                elif len(self.sequence[resid])<1:
                    inc -= 1
                elif len(self.sequence[resid])>1:
                    chain1 = re.split('(\d+)', resid)[0]
                    resnum1 = int(re.split('(\d+)', resid)[1])
                    present = False
                    for s in sym_sets:
                        if chain1+str(resnum1+inc) in s or (resnum1+inc)<1:
                            present = True
                            break
                    if not present:
                        new_set = {chain1+str(resnum1+inc)}
                        for elem in self.symmetric[resid]:
                            chain2 = re.split('(\d+)', elem)[0]
                            resnum2 = int(re.split('(\d+)', elem)[1])
                            new_set.add(chain2+str(resnum2+inc))
                        sym_sets.append(new_set)
                    
                        for aa in self.sequence[resid][1:]:
                            inc+=1
                            new_set = {chain1+str(resnum1+inc)}
                            for elem in self.symmetric[resid]:
                                chain2 = re.split('(\d+)', elem)[0]
                                resnum2 = int(re.split('(\d+)', elem)[1])
                                new_set.add(chain2+str(resnum2+inc))
                            sym_sets.append(new_set)
        new_sym = [list(s) for s in sym_sets]

        json_dict["symmetric"] = new_sym
        self.numbering = numbering
        self.jsondata = json_dict

    def _update_mutable_from_seq(self, seq):
        """updates mutable using newly generated sequence from mpnn"""
        i = 0
        newmut = {}
        seq = list(seq)
        for resid in self.sequence:
            if resid in self.mutable:
                new_res = copy.deepcopy(self.mutable[resid])
                if len(self.sequence[resid]) == 1:
                    #print(self.sequence[resid], new_res, seq, i)
                    new_res[0]["resid"][3] = seq[i]
                elif len(self.sequence[resid]) > 1:
                    #print(self.sequence[resid], new_res, seq, i)
                    for aa, j in zip(self.sequence[resid], range(len(self.sequence[resid]))):
                        if len(new_res)<=j:
                            new_res.append(copy.deepcopy(new_res[0]))
                            #print(new_res, j)
                            #print(new_res[j])
                            #print(new_res[j]["resid"])
                            new_res[j]["resid"][2] = j
                        new_res[j]["resid"][3] = aa

                newmut[resid] = new_res
            i+= len(self.sequence[resid])
            #print(i)
            if i>=len(seq):
                break
        return newmut

    def _update_sequence(self):
        """updates sequence dict based on changes to mutable"""
        #print(self.sequence)
        new_seq = {}
        for resid in self.sequence:
            if resid in self.mutable:
                pos_seq = ""
                for res in self.mutable[resid]:
                    pos_seq = pos_seq + res["resid"][3]

                if self.sequence[resid] != pos_seq:
                    self.sequence[resid] = pos_seq
        
    def _update_symmetric_positions(self, mutated):
        for mut_check in mutated:
            if mut_check in self.symmetric:
                for mut_id in self.symmetric[mut_check]:
                    chain = re.split('(\d+)', mut_id)[0]
                    index = int(re.split('(\d+)', mut_id)[1])
                    if len(self.mutable[mut_check]) != len(self.mutable[mut_id]):
                        self.mutable[mut_id] = copy.deepcopy(self.mutable[mut_check])
                        for pos in self.mutable[mut_id]:
                            pos["resid"][0] = chain
                            pos["resid"][1] = index

                    for res1, res2 in zip(self.mutable[mut_check], self.mutable[mut_id]):
                        if res1["resid"][3] != res2["resid"][3]:
                            #print("Before symmetric update", res1["resid"], res2["resid"])
                            res2["resid"][3] = res1["resid"][3]
                            #print("After symmetric update", res1["resid"], res2["resid"])
                        if res1["resid"][2] != res2["resid"][2]:
                            #print("Before symmetric update", res1["resid"], res2["resid"])
                            res2["resid"][2] = res1["resid"][2]
                            #print("After symmetric update", res1["resid"], res2["resid"])

    def _check_length_constraints(self):
        print("not working")

    def _check_symmetry(self):
        """checks symmetry"""
        for mut_id in self.mutable:
            if mut_id in self.symmetric:
                for sym_id in self.symmetric[mut_id]:
                    if self.sequence[mut_id] != self.sequence[sym_id]:
                        print(mut_id, sym_id, self.sequence[mut_id], self.sequence[sym_id])
                        print("not symmetric")

    def get_lengths(self, chains=None):
        lengths = []
        if chains:
            for chain in chains:
                lengths.append(len(self.jsondata["sequence"][chain]))

        else:
            for chain in self.jsondata["sequence"]:
                lengths.append(len(self.jsondata["sequence"][chain]))

        return lengths

    def get_sequence_string(self, divide=","):
        return divide.join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])

    def mutate(self, mut_percent = 0.125, num_mut_choice = [-1, 0, 1], var=0, var_weights = [0.8, 0.1, 0.1]):
        all_aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        #calculating number of mutants from mut_percent and maybe adding or subtracting one mutant for stochasticity
        #num_mut = round(len(str(self))*mut_percent) + random.choice([-1, 0, 1])
        num_mut=round(mut_percent*len(self.mutable.keys()))+random.choice(num_mut_choice)
        if num_mut<1:
            num_mut=1

        new_mut = copy.deepcopy(self.mutable) 

        mut_ids = random.sample(list(new_mut.keys()), len(new_mut.keys()))
        i = 0
        num_mut_curr = 0
        mutated = []
        while num_mut_curr < num_mut:
            mut_id = mut_ids[i]
            weights = new_mut[mut_id][0]["weights"]
            chain = re.split('(\d+)', mut_id)[0]
            aa = int(re.split('(\d+)', mut_id)[1])

            method = "sub"
            if var > 0:
                #print(len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])), len(self.sequence.keys()), var)
                if len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])) >= len(self.sequence.keys()) + var:
                    method = random.choices(["sub", "del"], [var_weights[0], var_weights[2]])[0]

                elif len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])) <= len(self.sequence.keys()) - var:
                    method = random.choices(["sub", "insert"], [var_weights[0], var_weights[1]])[0]

                else:
                    method = random.choices(["sub", "insert", "del"], var_weights)[0]

            print("mutating by", method, str(var), str(var_weights))
            old_res = new_mut[mut_id][-1]
            new_res = copy.deepcopy(old_res)
            if method == "sub":
                new_aa = random.choices(all_aas, weights)[0]
                if new_res["resid"][2] < 0:
                    print("Trying to mutate by substitution at a deletion. Mutating by insertion instead.")
                    new_res["resid"][2] = new_res["resid"][2] + 1

                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                # if substitution, replace dict with sub
                new_mut[mut_id][-1] = new_res
                mutated.append(mut_id)
            elif method == "insert":
                new_aa = random.choices(all_aas, weights)[0]
                new_res["resid"][2] = new_res["resid"][2] + 1
                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                #if insertion, append to list of dicts
                new_mut[mut_id].append(new_res)
                mutated.append(mut_id)
            elif method == "del":
                new_aa = ""
                new_res = copy.deepcopy(old_res)

                if new_res["resid"][2] < 0:
                    print("Trying to mutate by deletion at a deletion. Mutating by insertion instead.")
                    new_aa = random.choices(all_aas, weights)[0]
                    new_res["resid"][2] = new_res["resid"][2] + 1
                else:
                    new_res["resid"][2] = new_res["resid"][2] - 1

                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                # if deletion, replace dict with no-aa dict
                new_mut[mut_id][-1] = new_res
                mutated.append(mut_id)

            #print("after mutation", new_mut[mut_id])
            i+=1

        newseqobj = DesignSeq(sequence=self.sequence, mutable=new_mut, symmetric=self.symmetric)
        newseqobj._update_symmetric_positions(mutated)
        #newseqobj._check_symmetry()
        newseqobj._update_sequence()
        newseqobj._create_jsondata()
        newseqobj._check_symmetry()
        return newseqobj

    def _get_aa_identity(self, index):
        chain = re.split('(\d+)', index)[0]
        aa = int(re.split('(\d+)', index)[1])
        return self.seq[chain][aa-1]

    def crossover(self, otherDS, ncross = 1):
        chains = []
        for resid in self.mutable:
            chain = re.split('(\d+)', resid)[0]
            if chain not in chains:
                chains.append(chain)

        crossover_chain = random.choice(chains)
        mutated = []
        new_mut = {}
        
        mut_seq = [x for x in self.mutable if re.split('(\d+)', x)[0] == crossover_chain]
        other_mut_seq = [x for x in otherDS.mutable if re.split('(\d+)', x)[0] == crossover_chain]
        new_mut = copy.deepcopy(self.mutable)

        s=0
        points = random.sample(range(len(mut_seq)), ncross)
        for i, mut_id, other_mut_id in zip(range(len(mut_seq)), mut_seq, other_mut_seq):
            options = [self.mutable[mut_id], otherDS.mutable[other_mut_id]]
            if i in points:
                if s==0:
                    s=1
                else:
                    s=0
            new_mut[mut_id] = options[s]
            if s==1:
                mutated.append(mut_id)

        newseqobj = DesignSeq(sequence=self.sequence, mutable=new_mut, symmetric=self.symmetric)
        newseqobj._update_symmetric_positions(mutated)
        newseqobj._update_sequence()
        newseqobj._create_jsondata()
        newseqobj._check_symmetry()
        return newseqobj


class DesignSeqMSD(DesignSeq):
    def __init__(self, jsonfile=None, sequence=None, mutable=None, symmetric=[], jdata=None, seq=None):
        if jsonfile is not None:
            jseq, mut, sym, jdata = self._load_from_json(jsonfile)
            self.sequence = jseq
            self.mutable = mut
            self.symmetric = sym
            self.jsondata = jdata
            self._create_jsondata()

        elif seq is not None:
            self.sequence = sequence
            self.mutable = mutable
            self.mutable = self._update_mutable_from_seq(seq)
            self._update_sequence()
            self.symmetric = symmetric
            if jdata is not None:
                self.jsondata = jdata
            self._create_jsondata()
        else:
            self.sequence = sequence
            self.mutable = mutable
            self.symmetric = symmetric
            self._update_sequence()
            if jdata is not None:
                self.jsondata = jdata
            self._create_jsondata()
        
    def _create_jsondata(self):
        """creates json dictionary for dsobj that can be used as input to protein mpnn"""
        all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        new_des = []
        new_sym = []
        new_jsonseq = {}

        for resid in self.sequence:
            chain = re.split('(\d+)', resid)[0]
            if chain not in new_jsonseq:
                new_jsonseq[chain] = []

            new_jsonseq[chain].append(self.sequence[resid])

        newjsonseq_joined = {chain: "".join(new_jsonseq[chain]) for chain in new_jsonseq}
        json_dict = {"sequence":newjsonseq_joined}

        chains = []
        numbering = {}
        i = 1
        for resid in self.sequence:
            chain = re.split('(\d+)', resid)[0]
            #checking if we need to restart the numbering from 1 because new chain
            if chain not in chains:
                chains.append(chain)
                i = 1
                numbering[chain] = []

            resval = self.sequence[resid]
            if resid in self.mutable:
                for res in self.mutable[resid]:
                    if res["resid"][3] != '':
                        des_dict = {"chain":chain, "resid":i, "WTAA":res["resid"][3], "MutTo":res["MutTo"]}
                        new_des.append(des_dict)
                        #print("designable", des_dict)
                        numbering[chain].append(resid)
                        i+=1

                    else:
                        pass
            else:
                numbering[chain].append(resid) 
                i+=1


        json_dict["designable"] = new_des

        sym_sets = []
        chains = []
        inc = 0
        for resid in self.sequence:
            if resid in self.symmetric:
                chain1 = re.split('(\d+)', resid)[0]
                resnum1 = int(re.split('(\d+)', resid)[1])
                if len(self.sequence[resid]) == 1:
                    present = False
                    for s in sym_sets:
                        if chain1+str(resnum1+inc) in s or (resnum1+inc)<1:
                            present = True
                            break
                    if not present:
                        new_set = {chain1+str(resnum1+inc)}
                        for elem in self.symmetric[resid]:
                            chain2 = re.split('(\d+)', elem)[0]
                            resnum2 = int(re.split('(\d+)', elem)[1])
                            new_set.add(chain2+str(resnum2+inc))
                        sym_sets.append(new_set)
                elif len(self.sequence[resid])<1:
                    inc -= 1
                elif len(self.sequence[resid])>1:
                    chain1 = re.split('(\d+)', resid)[0]
                    resnum1 = int(re.split('(\d+)', resid)[1])
                    present = False
                    for s in sym_sets:
                        if chain1+str(resnum1+inc) in s or (resnum1+inc)<1:
                            present = True
                            break
                    if not present:
                        new_set = {chain1+str(resnum1+inc)}
                        for elem in self.symmetric[resid]:
                            chain2 = re.split('(\d+)', elem)[0]
                            resnum2 = int(re.split('(\d+)', elem)[1])
                            new_set.add(chain2+str(resnum2+inc))
                        sym_sets.append(new_set)
                    
                        for aa in self.sequence[resid][1:]:
                            inc+=1
                            new_set = {chain1+str(resnum1+inc)}
                            for elem in self.symmetric[resid]:
                                chain2 = re.split('(\d+)', elem)[0]
                                resnum2 = int(re.split('(\d+)', elem)[1])
                                new_set.add(chain2+str(resnum2+inc))
                            sym_sets.append(new_set)
        new_sym = [list(s) for s in sym_sets]

        json_dict["symmetric"] = new_sym
        self.numbering = numbering
        if self.jsondata is not None:
            if "tied_betas" in self.jsondata:
                json_dict["tied_betas"] = self.jsondata["tied_betas"]
            if "chain_key" in self.jsondata:
                json_dict["chain_key"] = self.jsondata["chain_key"]
        json_dict["chain_key"] = self.jsondata["chain_key"]
        self.jsondata = json_dict

    def crossover(self, otherDS, ncross = 1):
        chains = []
        for resid in self.mutable:
            chain = re.split('(\d+)', resid)[0]
            if chain not in chains:
                chains.append(chain)

        crossover_chain = random.choice(chains)
        mutated = []
        new_mut = {}
        
        mut_seq = [x for x in self.mutable if re.split('(\d+)', x)[0] == crossover_chain]
        other_mut_seq = [x for x in otherDS.mutable if re.split('(\d+)', x)[0] == crossover_chain]
        new_mut = copy.deepcopy(self.mutable)

        s=0
        points = random.sample(range(len(mut_seq)), ncross)
        for i, mut_id, other_mut_id in zip(range(len(mut_seq)), mut_seq, other_mut_seq):
            options = [self.mutable[mut_id], otherDS.mutable[other_mut_id]]
            if i in points:
                if s==0:
                    s=1
                else:
                    s=0
            new_mut[mut_id] = options[s]
            if s==1:
                mutated.append(mut_id)

        newseqobj = DesignSeqMSD(sequence=self.sequence, mutable=new_mut, symmetric=self.symmetric, jdata=self.jsondata)
        newseqobj._update_symmetric_positions(mutated)
        newseqobj._update_sequence()
        newseqobj._create_jsondata()
        #newseqobj._check_symmetry()
        return newseqobj

    def mutate(self, mut_percent = 0.125, num_mut_choice = [-1, 0, 1], var=0, var_weights = [0.8, 0.1, 0.1]):
        all_aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        #calculating number of mutants from mut_percent and maybe adding or subtracting one mutant for stochasticity
        #num_mut = round(len(str(self))*mut_percent) + random.choice([-1, 0, 1])
        num_mut=round(mut_percent*len(self.mutable.keys()))+random.choice(num_mut_choice)
        if num_mut<1:
            num_mut=1

        new_mut = copy.deepcopy(self.mutable) 

        mut_ids = random.sample(list(new_mut.keys()), len(new_mut.keys()))
        i = 0
        num_mut_curr = 0
        mutated = []
        while num_mut_curr < num_mut:
            mut_id = mut_ids[i]
            weights = new_mut[mut_id][0]["weights"]
            chain = re.split('(\d+)', mut_id)[0]
            aa = int(re.split('(\d+)', mut_id)[1])

            method = "sub"
            if var > 0:
                #print(len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])), len(self.sequence.keys()), var)
                if len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])) >= len(self.sequence.keys()) + var:
                    method = random.choices(["sub", "del"], [var_weights[0], var_weights[2]])[0]

                elif len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])) <= len(self.sequence.keys()) - var:
                    method = random.choices(["sub", "insert"], [var_weights[0], var_weights[1]])[0]

                else:
                    method = random.choices(["sub", "insert", "del"], var_weights)[0]

            old_res = new_mut[mut_id][-1]
            new_res = copy.deepcopy(old_res)
            if method == "sub":
                new_aa = random.choices(all_aas, weights)[0]
                if new_res["resid"][2] < 0:
                    print("Trying to mutate by substitution at a deletion. Mutating by insertion instead.")
                    new_res["resid"][2] = new_res["resid"][2] + 1

                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                # if substitution, replace dict with sub
                new_mut[mut_id][-1] = new_res
                mutated.append(mut_id)
            elif method == "insert":
                new_aa = random.choices(all_aas, weights)[0]
                new_res["resid"][2] = new_res["resid"][2] + 1
                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                #if insertion, append to list of dicts
                new_mut[mut_id].append(new_res)
                mutated.append(mut_id)
            elif method == "del":
                new_aa = ""
                new_res = copy.deepcopy(old_res)

                if new_res["resid"][2] < 0:
                    print("Trying to mutate by deletion at a deletion. Mutating by insertion instead.")
                    new_aa = random.choices(all_aas, weights)[0]
                    new_res["resid"][2] = new_res["resid"][2] + 1
                else:
                    new_res["resid"][2] = new_res["resid"][2] - 1

                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                # if deletion, replace dict with no-aa dict
                new_mut[mut_id][-1] = new_res
                mutated.append(mut_id)

            #print("after mutation", new_mut[mut_id])
            i+=1

        newseqobj = DesignSeqMSD(sequence=self.sequence, mutable=new_mut, symmetric=self.symmetric, jdata=self.jsondata)
        newseqobj._update_symmetric_positions(mutated)
        #newseqobj._check_symmetry()
        newseqobj._update_sequence()
        newseqobj._create_jsondata()
        return newseqobj


if __name__ == "__main__":
    dsobj = DesignSeq(jsonfile="/pine/scr/a/m/amritan/kuhlmanlab/run_evopro/tests/test01/residue_specs.json")
    print(dsobj.mutable)
    dupdsobj = copy.deepcopy(dsobj)
    seq = "DLLRRMLGMVIRMLGVFTKLLGKILMIPAGIYAPICVTVRYFETVGEALERAGILLRGRDRAGKPRLTPAAREILKEALKAAEEAVDVLTLDITNITTSHQRKMESLNFIRAHTPYINIYNCEPANPSEKNSPLMQYCKALQNLRLAVLNVGLEIAKLAVKLISADLLRRMLGMVIRMLGVFTKLLGKILMIPAGIYAPICVTVRYFETVGEALERAGILLRGRDRAGKPRLTPAAREILKEALKAAEEAVDVLTLDITNITTSHQRKMESLNFIRAHTPYINIYNCEPANPSEKNSPLMQYCKALQNLRLAVLNVGLEIAKLAVKLISA"
    newdsobj = DesignSeq(seq=seq, sequence=dsobj.sequence, mutable=dsobj.mutable, symmetric=dsobj.symmetric)
    print(newdsobj.mutable)
    print(dsobj == newdsobj)
    print(dsobj == dupdsobj)
    """
    print(dsobj.sequence)
    newdsobj = dsobj.mutate(var=5, var_weights = [0.1, 0.8, 0.1])
    print(newdsobj.sequence)
    #print(dsobj.sequence, newdsobj.sequence)
    #print("after", newdsobj.sequence, newdsobj.mutable, newdsobj.symmetric,newdsobj.jsondata)
    #print(newdsobj.sequence, newdsobj.jsondata["sequence"])
    newdsobj2 = newdsobj.mutate(var=5, var_weights = [0, 0.8, 0.1])
    print(newdsobj2.sequence)
    newdsobj3 = newdsobj2.mutate(var=5, var_weights = [0, 0.8, 0.1])
    print(newdsobj3.sequence)
    print(dsobj.jsondata["sequence"], newdsobj.jsondata["sequence"], newdsobj2.jsondata["sequence"], newdsobj3.jsondata["sequence"])
    #print(dsobj.get_lengths())
    """
