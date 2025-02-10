from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.cealign import CEAligner

import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.score_funcs.calculate_rmsd import kabsch_rmsd, kabsch_rmsd_superimposeall
from evopro.utils.pdb_parser import get_coordinates_pdb

import math
import pickle
import numpy as np
import re

def calculate_rmsd_cealigner(pdb_file_1, pdb_file_2):
    # Initialize the PDB parser
    parser = PDBParser()

    # Parse the two structures
    structure1 = parser.get_structure('structure1', pdb_file_1)
    structure2 = parser.get_structure('structure2', pdb_file_2)

    # Initialize the CEAligner
    ce_aligner = CEAligner()

    # Perform the alignment
    ce_aligner.set_reference(structure1)
    ce_aligner.align(structure2)

    # Get the RMSD from the alignment
    rmsd = ce_aligner.rms
    
    io = PDBIO()
    io.set_structure(structure1)
    io.save("/work/users/a/m/amritan/evopro_tests/rmsd/for_nikka/out_1.pdb")
    io.set_structure(structure2)
    io.save("/work/users/a/m/amritan/evopro_tests/rmsd/for_nikka/out_2.pdb")

    return rmsd

def calculate_rmsd_evopro(pdbfile1, pdbfile2, ca_only=True, translate=True):
    with open(pdbfile1, "r") as f:
        pdb1_string = f.read()
        
    with open(pdbfile2, "r") as f:
        pdb2_string = f.read()

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1_string)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2_string)
    reslist1 = [x for x in residues1.keys()]
    reslist2 = [x for x in residues2.keys()]
    
    A = []
    for res in reslist1:
        for atom in residues1[res]:
            if ca_only:
                if atom[1] == 'CA':
                    A.append(list(atom[-1]))
            else:
                A.append(list(atom[-1]))
    B = []
    for res in reslist2:
        for atom in residues2[res]:
            if ca_only:
                if atom[1] == 'CA':
                    B.append(list(atom[-1]))
            else:
                B.append(list(atom[-1]))
    A = np.array(A)
    B = np.array(B)
    A = A.astype(float)
    B = B.astype(float)
    rmsd = kabsch_rmsd(A, B, translate=translate)
    return rmsd

if __name__ == "__main__":
    
    print(calculate_rmsd_cealigner("/work/users/a/m/amritan/evopro_tests/rmsd/for_nikka/design_56.pdb", "/work/users/a/m/amritan/evopro_tests/rmsd/for_nikka/seq_0_final_model_1_chainAB.pdb"))    
    print(calculate_rmsd_evopro("/work/users/a/m/amritan/evopro_tests/rmsd/for_nikka/design_56.pdb", "/work/users/a/m/amritan/evopro_tests/rmsd/for_nikka/seq_0_final_model_1_chainAB.pdb"))