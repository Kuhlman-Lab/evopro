import sys, importlib
from scoring.standard_score_funcs import *
from utils.parsing_utils import *
from utils.parsing import parse_results
from utils.hbonding import score_hydrogen_bonds
from utils.utils import check_omegaconf_key_exists

def score_overall(results, sequences, dsobj, conf):
    all_chains = [x.strip() for x in conf.structure_prediction.structure_pred_chains]

    label = ""
    num_chains_seq = len(str(dsobj).split(","))
    for i in range(num_chains_seq-1):
        label += "," 
    results, pdbs = parse_results(results, conf)
    label += ",overall_score"
    plddt_scores = []
    if conf.scoring.plddt.score_plddt:
        plddt_chains = [x.strip() for x in conf.scoring.plddt.preds.split(",")]
        plddt_weights = parse_weights(conf.scoring.plddt.weights, plddt_chains)
        residues = get_residues(conf.scoring.plddt.residues, dsobj)

        i = 0
        for result, pred in zip(results, all_chains):
            if pred in plddt_chains:
                if 'plddt' in result:
                    plddt_score = get_avg_plddt(result, reslist=residues)
                else:
                    residue_plddt_df, plddt_score = calculate_residue_plddt(result, reslist=residues)
                plddt_scores.append((plddt_score, plddt_weights[i]))
                i+=1
                label += ",pLDDT_chain" + pred + ",weighted_pLDDT_chain" + pred
                
    ptm_scores = []
    if conf.scoring.ptm.score_ptm:
        ptm_chains = [x.strip() for x in conf.scoring.ptm.preds.split(",")]
        ptm_weights = parse_weights(conf.scoring.ptm.weights, ptm_chains)
        i=0
        for result, pred in zip(results, all_chains):
            if pred in ptm_chains:
                if 'ptm' in result:
                    ptm = result['ptm']
                    if type(ptm) != float:
                        ptm = float(ptm)
                    ptm_scores.append((ptm, ptm_weights[i]))
                    i+=1
        label += ",pTM,weighted_pTM"

    if conf.scoring.iptm.score_iptm:
        iptm_chains = [x.strip() for x in conf.scoring.iptm.preds.split(",")]
        iptm_weights = parse_weights(conf.scoring.iptm.weights, iptm_chains)
        i=0
        for result, pred in zip(results, all_chains):
            if pred in iptm_chains:
                if 'iptm' in result:
                    iptm = result['iptm']
                    if type(iptm) != float:
                        iptm = float(iptm)
                    ptm_scores.append((iptm, iptm_weights[i]))
                    i+=1
        label += ",ipTM,weighted_ipTM"
                
    contacts_scores = []
    if conf.scoring.contacts.score_contacts:
        contacts_chains = [x.strip() for x in conf.scoring.contacts.preds.split(",")]
        contacts_weights = parse_weights(conf.scoring.contacts.weights, contacts_chains)
        reslist1 = get_tokens(conf.scoring.contacts.interface_1, dsobj)
        reslist2 = get_tokens(conf.scoring.contacts.interface_2, dsobj)
        label += ",ContactScore,weighted_ContactScore"
        
        if conf.scoring.contacts.bonus_residues:
            bonus_residues = get_tokens(conf.scoring.contacts.bonus_residues, dsobj)
            label += ",BonusScore,weighted_BonusScore"
        if conf.scoring.contacts.penalty_residues:
            penalty_residues = get_tokens(conf.scoring.contacts.penalty_residues, dsobj)
            label += ",PenaltyScore,weighted_PenaltyScore"
        
        i=0
        for result, seq, pred in zip(results, sequences, all_chains):
            if pred in contacts_chains:
                
                pairs, contacts_score = score_contacts_pae_weighted(result, seq, reslist1, reslist2, conf)
                contacts_scores.append((contacts_score, contacts_weights[i]))
                i+=1
                
                if conf.scoring.contacts.bonus_residues:
                    bonus_contact, bonus = score_contact(result, bonus_residues, reslist2, dist=conf.scoring.contacts.bonus_distance)
                    contacts_scores.append((bonus, conf.scoring.contacts.bonus_weight))
                if conf.scoring.contacts.penalty_residues:
                    penalty_contact, penalty = score_contact(result, penalty_residues, reslist2, dist=conf.scoring.contacts.penalty_distance)
                    contacts_scores.append((penalty, conf.scoring.contacts.penalty_weight))
                    
    rmsd_scores = []
    #conformational difference score 
    if conf.scoring.conf_change.score_conf_change:
        conf_chains = [x.strip() for x in conf.scoring.conf_change.preds.split(",")]
        pred1 = conf_chains[0]
        pred2 = conf_chains[1]
        conf_change_weight = conf.scoring.conf_change.weights
        for result, pred in zip(results, all_chains):
            if pred == pred1:
                pdb1 = result['pdb']
            if pred == pred2:
                pdb2 = result['pdb']
        
        reslist = get_residues(conf.scoring.conf_change.residues, dsobj)
        conf_change = score_rmsd(pdb1, pdb2, reslist=reslist)
        rmsd_scores.append((conf_change, conf_change_weight))
        label += ",ConfDiffScore,weighted_ConfDiffScore"
        
        if conf.scoring.conf_change.threshold_penalty:
            if conf_change > conf.scoring.conf_change.threshold:
                penalty = calc_penalty_score(conf_change, conf.scoring.conf_change.threshold, conf.scoring.conf_change.spring_constant)
                rmsd_scores.append((conf_change - conf.scoring.conf_change.threshold, penalty))
                label += ",ConfDiffScoreAboveThreshold,penalty_ConfDiffScoreAboveThreshold"
        
    #rmsd score
    if conf.scoring.rmsd.score_rmsd:
        rmsd_chains = [x.strip() for x in conf.scoring.rmsd.preds.split(",")]
        rmsd_weights = parse_weights(conf.scoring.rmsd.weights, rmsd_chains)
        reslist = get_residues(conf.scoring.rmsd.residues, dsobj)
        
        i=0
        for result, pred in zip(results, all_chains):
            if pred in rmsd_chains:
                with open(conf.scoring.rmsd.rmsd_pdb, "r") as f:
                    rmsd_pdb = f.read()
                rmsd_score = score_rmsd(result['pdb'], rmsd_pdb, reslist)
                rmsd_scores.append((rmsd_score, rmsd_weights[i]))
                i+=1
        label += ",RMSDScore,weighted_RMSDScore"
        if conf.scoring.rmsd.threshold_penalty:
            for rmsd_score in rmsd_scores:
                if rmsd_score[0] > conf.scoring.rmsd.threshold:
                    penalty = calc_penalty_score(rmsd_score[0], conf.scoring.rmsd.threshold, conf.scoring.rmsd.spring_constant)
                    rmsd_scores.append((rmsd_score[0] - conf.scoring.rmsd.threshold, penalty))
                    label += ",RMSDScoreAboveThreshold,penalty_RMSDScoreAboveThreshold"
        
    other_scores = []
    if check_omegaconf_key_exists("scoring.hbond", conf):
        if conf.scoring.hbond.score_hbond:
            hbond_chains = [x.strip() for x in conf.scoring.hbond.preds.split(",")]
            hbond_weights = parse_weights(conf.scoring.hbond.weights, hbond_chains)
            reslist1 = get_tokens(conf.scoring.contacts.interface_1, dsobj)
            reslist2 = get_tokens(conf.scoring.contacts.interface_2, dsobj)

            i=0
            for result, pred in zip(results, all_chains):
                if pred in hbond_chains:
                    pairs, hbond_score = score_hydrogen_bonds(result, reslist1, reslist2, max_don_h_dist=conf.scoring.hbond.max_don_h_dist)
                    other_scores.append((hbond_score, hbond_weights[i]))
                    i+=1
                    print("Hbond pairs:", pairs)

            label += ",HBondScore,weighted_HBondScore"
    
    if check_omegaconf_key_exists("scoring.hydrophobic_interface", conf):
        if conf.scoring.hydrophobic_interface.score_hydrophobic_interface:
            hydphob_chains = [x.strip() for x in conf.scoring.hydrophobic_interface.preds.split(",")]
            hydphob_weights = parse_weights(conf.scoring.hydrophobic_interface.weights, hydphob_chains)
            t = "interface"
            i=0
            for result, pred in zip(results, all_chains):
                if pred in hydphob_chains:
                    hydphob_score = score_hydrophobicity(result, conf, score_type=t)
                    print("Hydrophobicity score:", hydphob_score)
                    other_scores.append((hydphob_score, hydphob_weights[i]))
                    i+=1
            label += ",HydrophobicityScore,weighted_HydrophobicityScore"
            
    if check_omegaconf_key_exists("scoring.hydrophobic_surface", conf):
        if conf.scoring.hydrophobic_surface.score_hydrophobic_surface:
            hydphob_chains = [x.strip() for x in conf.scoring.hydrophobic_surface.preds.split(",")]
            hydphob_weights = parse_weights(conf.scoring.hydrophobic_surface.weights, hydphob_chains)
            t = "surface"
            i=0
            for result, pred in zip(results, all_chains):
                if pred in hydphob_chains:
                    hydphob_score = score_hydrophobicity(result, conf, score_type=t)
                    print("Hydrophobicity score:", hydphob_score)
                    other_scores.append((hydphob_score, hydphob_weights[i]))
                    i+=1
            label += ",HydrophobicityScore,weighted_HydrophobicityScore"
            
    custom_scores = []
    if conf.scoring.custom.score_custom:
        custom_chains = [x.strip() for x in conf.scoring.custom.preds.split(",")]
        custom_weights = parse_weights(conf.scoring.custom.weights, custom_chains)
        i=0
        for result, pred in zip(results, all_chains):
            if pred in custom_chains:
                custom_score = score_custom(result, conf, dsobj)
                custom_scores.append((custom_score, custom_weights[i]))
                i+=1
        label += ",CustomScore,weighted_CustomScore"
    
    overall_score = 0
    all_scores = plddt_scores + ptm_scores + contacts_scores + rmsd_scores + other_scores + custom_scores
    print(all_scores)
    for s in all_scores:
        overall_score += s[0]*s[1]
    
    all_scores.insert(0, overall_score)
    return all_scores, pdbs, results, label

def score_custom(result, conf, dsobj):
    try:
        if conf.scoring.custom:
            scorefile = conf.scoring.custom.custom_script.rsplit("/", 1)
            scorepath = scorefile[0]
            scorefilename = scorefile[1].split(".")[0]

            sys.path.append(scorepath)
            mod = importlib.import_module(scorefilename)
            scorefunc = getattr(mod, conf.scoring.custom.custom_score_func)
            return scorefunc(result, conf, dsobj)
    except:
        raise ValueError("Invalid score function")
    
def calc_penalty_score(score, threshold, spring_constant):
    """
    Calculate the penalty score based on the threshold and spring constant.
    """
    if score > threshold:
        penalty = spring_constant * (score - threshold)
        return penalty
    else:
        return 0