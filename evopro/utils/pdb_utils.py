from objects.chemical import to1letter  
from utils.parsing_utils import get_coordinates_pdb

def get_seq_from_pdb(pdb_str):
    seqs = []
    chains = []

    seq = []
    pdb_split = pdb_str.split("\n")
    pdb_split = [x for x in pdb_split if x]
    for lin in pdb_split:
        l = lin.strip().split()
        if 'ATOM' in l[0] or 'HETATM' in l[0]:
            if l[2] == "CA":
                if l[4] not in chains:
                    chains.append(l[4])
                    seqs.append("".join(seq))
                    seq = []
                
                seq.append(to1letter[l[3]])
    seqs.append("".join(seq))
    
    seqs = [x for x in seqs if x]
    return seqs

def change_chainid_pdb(pdb, old_chain="A", new_chain="B"):
    pdb_lines = [x for x in pdb.split("\n") if x]
    #print(pdb_lines)
    new_pdb_lines = []
    for lin in pdb_lines:
        lin_list = list(lin.strip())
        if len(lin_list)>21:
            if lin_list[21] == old_chain:
                lin_list[21] = new_chain
        new_lin = "".join(lin_list) + "\n"
        new_pdb_lines.append(new_lin)
    #print(pdb_lines, new_pdb_lines)
    return "".join(new_pdb_lines)

def transform_pdb_location(pdb, offset_vals = (0,0,0)):
    pdb_lines = [x for x in pdb.split("\n") if x]
    #print(pdb_lines)
    new_pdb_lines = []
    for lin in pdb_lines:
        lin_list = list(lin.strip())
        if lin.startswith("ATOM") or lin.startswith("HETATM"):
            x = lin[30:38].strip(' ')
            y = lin[38:46].strip(' ')
            z = lin[46:54].strip(' ')
            x = str(round(float(x) + offset_vals[0], 2))
            y = str(round(float(y) + offset_vals[1], 2))
            z = str(round(float(z) + offset_vals[2], 2))
            while len(x)<8:
                x = " " + x
            while len(y)<8:
                y = " " + y
            while len(z)<8:
                z = " " + z
            lin_list[30:38] = list(x)
            lin_list[38:46] = list(y)
            lin_list[46:54] = list(z)
            
        new_lin = "".join(lin_list) + "\n"
        new_pdb_lines.append(new_lin)
    #print(pdb_lines, new_pdb_lines)
    return "".join(new_pdb_lines)


def find_max_coordinates(pdb):
    chains, residues, resindices = get_coordinates_pdb(pdb)
    all_atoms = []
    all_x = []
    all_y = []
    all_z = []
    for res in residues:
        #print(res, residues[res])
        for atom in residues[res]:
            #print(atom)
            all_atoms.append(atom[-1])
            all_x.append(float(atom[-1][0]))
            all_y.append(float(atom[-1][1]))
            all_z.append(float(atom[-1][2]))
    
    return (max(all_x), max(all_y), max(all_z))
    #print(residues)

def find_min_coordinates(pdb):
    chains, residues, resindices = get_coordinates_pdb(pdb)
    all_atoms = []
    all_x = []
    all_y = []
    all_z = []
    for res in residues:
        #print(res, residues[res])
        for atom in residues[res]:
            #print(atom)
            all_atoms.append(atom[-1])
            all_x.append(float(atom[-1][0]))
            all_y.append(float(atom[-1][1]))
            all_z.append(float(atom[-1][2]))
    
    return (min(all_x), min(all_y), min(all_z))

def append_pdbs(pdb1, pdb2):
    # with open("pdb1.pdb", "w") as pdb1f:
    #     pdb1f.write(pdb1)
    # with open("pdb2.pdb", "w") as pdb2f:
    #     pdb2f.write(pdb2)
    pdb1_lines = [x.strip() for x in pdb1.split("\n") if x]
    new_pdb_lines = []
    for lin in pdb1_lines:
        if lin.startswith("END"):
            break
        new_pdb_lines.append(lin+"\n")
        
    max_pdb1 = find_max_coordinates(pdb1)
    min_pdb2 = find_min_coordinates(pdb2)
    
    offset_vals = [abs(a) + abs(b) + 100 for a,b in zip(max_pdb1, min_pdb2)]
    #print(tuple(offset_vals))
    transformed_pdb2 = transform_pdb_location(pdb2, tuple(offset_vals))
    
    pdb2_lines = [x.strip() for x in transformed_pdb2.split("\n") if x]
    numbering = int(new_pdb_lines[-1][6:11].strip()) + 1
    for lin in pdb2_lines:
        lin_list = list(lin.strip())
        number = str(numbering)
        if lin.startswith("ATOM") or lin.startswith("HETATM") or lin.startswith("TER"):
            while len(number)<5:
                number = " " + number
            lin_list[6:11] = list(number)
            numbering+=1
        new_lin = "".join(lin_list) + "\n"
        if not new_lin.startswith("MODEL"):
            new_pdb_lines.append(new_lin)
    return "".join(new_pdb_lines)