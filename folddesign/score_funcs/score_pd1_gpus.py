from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.score_funcs.score_funcs_of import score_conf_of, score_confidence_pairs_of, score_confidence_interface_of
from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
import os
import subprocess
import shutil

pd1 = "/proj/kuhl_lab/folddesign/folddesign/data/pd1.pdb"

def score_pd1_of(results, pdb):
    pocketresidues = [('31','SER'),('32','PHE'),('33','HIS'),('34','VAL'),('35','VAL'),('36','TRP'),('37','HIS'),('38','ARG'),('39','GLU'),('40','SER'),('41','PRO'),('42','SER'),('43','GLY'),('44','GLN'),('45','THR'),('46','ASP'),('47','THR'),('48','LEU'),('49','ALA'),('50','ALA'),('51','PHE'),('52','PRO'),('53','GLU'),('54','ASP'),('55','ARG'),('56','SER'),('57','GLN'),('58','PRO'),('88','GLY'),('89','THR'),('90','TYR'),('91','VAL'),('92','CYS'),('93','GLY'),('94','VAL'),('95','ILE'),('96','SER'),('102','GLN'),('103','ILE'),('104','LYS'),('105','GLU'),('106','SER'),('107','LEU'),('108','ARG')]
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False, of=True)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>120]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, of=True)
    confidencescore = score_confidence_interface_of(results, reslist1, reslist2, resindices)
    confscore = score_conf_of(results, reslist2, resindices)
    score = -confscore + confidencescore
    return score, (confidencescore, confscore), pdb, results

def score_pd1_func(results, orient=None):
    from alphafold.common import protein
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
    from alphafold.common import protein
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
    from alphafold.common import protein
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
    from alphafold.common import protein
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
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    score = -confscore2
    return score, (confscore2), pdb, results

def score_binder_rmsd(pdb1, pdb2, with_linker=False):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1, fil = False)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2, fil = False)
    if len(chains1)>1:
        reslist1 = [x for x in residues1.keys() if x.startswith("B")]
    else:
        reslist1 = [x for x in residues1.keys() if int(x.split("_")[2])>135]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2)

    rmsd = rmsd_binder
    if with_linker:
        with open(pd1, "r") as pd1f:
            pd1_pdb = pd1f.read()
        chains3, residues3, resindices3 = get_coordinates_pdb(pd1_pdb, fil = False)
        reslist3 = [x for x in residues3.keys()]
        reslist4 = [x for x in residues1.keys() if int(x.split("_")[2])<120]
        rmsd_pd1 = get_rmsd(reslist3, pd1_pdb, reslist4, pdb1)
        rmsd = rmsd + rmsd_pd1
    return rmsd


def score_pd1_func_6(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(chains)>1:
        return score_pd1_func_3(results)
    else:
        return score_pd1_func_5(results)

def score_pd1_func_7(results):
    #with orientation restraint
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(chains)>1:
        return score_pd1_func_4(results)
    else:
        return score_pd1_func_5(results)

def score_pd1_func_8(results):
    #with linker
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(residues.keys())>135:
        return score_pd1_func_2(results)
    else:
        return score_pd1_func_5(results)

def score_pd1_otherinterface(results):
    from alphafold.common import protein
    pocketresidues = [('5', 'THR'), ('7', 'SER'), ('8', 'PRO'), ('11', 'LEU'), ('17', 'ASP'), ('18', 'ASN'), ('20', 'THR'), ('22', 'THR'), ('24', 'SER'), ('65', 'ARG'), ('67', 'THR'), ('70', 'PRO'), ('71', 'ASN'), ('73', 'ARG'), ('74', 'ASP'), ('76', 'HIS'), ('78', 'SER'), ('80', 'VAL'), ('81', 'ARG')]
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = ["A_"+x[1]+"_"+x[0] for x in pocketresidues]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdb, results

def score_pd1_func_9(results):
    #trying to design binder on some other interface (not the pocket) as a negative control
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(chains)>1:
        return score_pd1_otherinterface(results)
    else:
        return score_pd1_func_5(results)

if __name__ == "__main__":
    results = []
    score_pd1_func(results)
