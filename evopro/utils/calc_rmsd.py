import copy
import numpy as np
from evopro.utils.write_pdb import PDBio

def fit_rms(ref_c,c):
    # move geometric center to the origin
    ref_trans = np.average(ref_c, axis=0)
    ref_c = ref_c - ref_trans
    c_trans = np.average(c, axis=0)
    c = c - c_trans

    # covariance matrix
    C = np.dot(c.T, ref_c)

    # Singular Value Decomposition
    (r1, s, r2) = np.linalg.svd(C)

    # compute sign (remove mirroring)
    if np.linalg.det(C) < 0:
        r2[2,:] *= -1.0
    U = np.dot(r1, r2)
    return (c_trans, U, ref_trans)

class RMSDcalculator:
    def __init__(self, atoms1, atoms2, name=None):
        xyz1 = self.get_xyz(atoms1, name=name)
        xyz2 = self.get_xyz(atoms2, name=name)
        self.set_rmsd(xyz1, xyz2)

    def get_xyz(self, atoms, name=None):
        xyz = []
        for atom in atoms:
            if name:
                if atom.name != name: continue
            xyz.append([atom.x, atom.y, atom.z])
        return np.array(xyz)

    def set_rmsd(self, c1, c2):
        self.rmsd = 0.0
        self.c_trans, self.U, self.ref_trans = fit_rms(c1, c2)
        new_c2 = np.dot(c2 - self.c_trans, self.U) + self.ref_trans
        self.rmsd = np.sqrt( np.average( np.sum( ( c1 - new_c2 )**2, axis=1 ) ) )

    def get_aligned_coord(self, atoms, name=None):
        new_c2 = copy.deepcopy(atoms)
        for atom in new_c2:
            atom.x, atom.y, atom.z = np.dot(np.array([atom.x, atom.y, atom.z]) - self.c_trans, self.U) + self.ref_trans
        return new_c2

if __name__ == '__main__':
    pdbf1 = '../tests/pd1_threehelix_run1/sequence_0_model_1_unrelaxed.pdb'; pdbf2 = '../tests/pd1_threehelix_run1/sequence_10_model_1_unrelaxed.pdb'
    with open(pdbf1, "r") as f1:
        pdbs1 = f1.read()
    with open(pdbf2, "r") as f2:
        pdbs2 = f2.read()
    pdb1 = PDBio(pdbs1); pdb2 = PDBio(pdbs2)
    atoms1 = pdb1.get_atoms(to_dict=False); atoms2 = pdb2.get_atoms(to_dict=False)

    RMSDcalculator = RMSDcalculator(atoms1, atoms2, name='CA')
    rmsd = RMSDcalculator.rmsd
    new_atoms = RMSDcalculator.get_aligned_coord(atoms2)
    #pdb2.write_pdb('aligned_%s' % pdbf2, new_atoms)
    print('RMSD : %8.3f' % rmsd)
    #print('New structure file: ', 'aligned_%s' % pdbf2)
