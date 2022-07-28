from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
from folddesign.utils.pdb_parser import get_coordinates_pdb
import os
import subprocess
from alphafold.common import protein
import shutil

def score_pd1_func(results, orient=None):
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False, orient=orient)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    score = -contactscore*100 + confidencescore/10.0 + confscore2
    return score, (contactscore, confidencescore, confscore2), contacts, pdb, results

def score_pd1_func_2(results, orient=None):
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('97','LEU'),('98','ALA'),('99','PRO'),('100','LYS'),('101','ILE'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>135]
    #reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False, orient=None)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdb, results

def score_pd1_func_3(results):
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('97','LEU'),('98','ALA'),('99','PRO'),('100','LYS'),('101','ILE'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdb, results

def score_pd1_func_4(results):
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('97','LEU'),('98','ALA'),('99','PRO'),('100','LYS'),('101','ILE'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    orient = [reslist2[0], ['A_GLY_43']]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False, orient=orient)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdb, results

def score_pd1_func_5(results, orient=None):
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    score = -confscore2
    return score, (confscore2), pdb, results

def score_binder_rmsd(pdb1, pdb2):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1, fil = False)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2, fil = False)
    if len(chains1)>1:
        reslist1 = [x for x in residues1.keys() if x.startswith("B")]
    else:
        reslist1 = [x for x in residues1.keys() if int(x.split("_")[2])>135]
    reslist2 = [x for x in residues2.keys()]
    rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2)
    return rmsd


def score_pd1_func_6(results):
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(chains)>1:
        return score_pd1_func_3(results)
    else:
        return score_pd1_func_5(results)

def score_pd1_func_7(results):
    #with orientation restraint
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(chains)>1:
        return score_pd1_func_4(results)
    else:
        return score_pd1_func_5(results)

def score_pd1_func_8(results):
    #with linker
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(residues.keys())>135:
        return score_pd1_func_2(results)
    else:
        return score_pd1_func_5(results)

if __name__ == "__main__":
    results = []
    score_pd1_func(results)
