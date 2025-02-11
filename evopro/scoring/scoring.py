import sys, importlib
from scoring.standard_score_funcs import *
from utils.parsing_utils import *
from utils.parsing import parse_results

def score_overall(results, sequences, dsobj, conf):
    all_chains = [x.strip() for x in conf.structure_prediction.structure_pred_chains]

    results, pdbs = parse_results(results, conf)

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
                
    ptm_scores = []
    if conf.scoring.ptm.score_ptm:
        ptm_chains = [x.strip() for x in conf.scoring.ptm.preds.split(",")]
        ptm_weights = parse_weights(conf.scoring.ptm.weights, ptm_chains)
        i=0
        for result, pred in zip(results, all_chains):
            if pred in ptm_chains:
                if 'ptm' in result:
                    ptm_scores.append((result['ptm'], ptm_weights[i]))
                    i+=1
                    
    if conf.scoring.iptm.score_iptm:
        iptm_chains = [x.strip() for x in conf.scoring.iptm.preds.split(",")]
        iptm_weights = parse_weights(conf.scoring.iptm.weights, iptm_chains)
        i=0
        for result, pred in zip(results, all_chains):
            if pred in iptm_chains:
                if 'iptm' in result:
                    ptm_scores.append((result['iptm'], iptm_weights[i]))
                    i+=1
                
    contacts_scores = []
    if conf.scoring.contacts.score_contacts:
        contacts_chains = [x.strip() for x in conf.scoring.contacts.preds.split(",")]
        contacts_weights = parse_weights(conf.scoring.contacts.weights, contacts_chains)
        reslist1 = get_tokens(conf.scoring.contacts.interface_1, dsobj)
        reslist2 = get_tokens(conf.scoring.contacts.interface_2, dsobj)
        
        if conf.scoring.contacts.bonus_residues:
            bonus_residues = get_tokens(conf.scoring.contacts.bonus_residues, dsobj)
        if conf.scoring.contacts.penalty_residues:
            penalty_residues = get_tokens(conf.scoring.contacts.penalty_residues, dsobj)
        
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
            # print(len(result))
            # print(result.keys())
            if pred == pred1:
                pdb1 = result['pdb']
            if pred == pred2:
                pdb2 = result['pdb']
        
        reslist = get_residues(conf.scoring.conf_change.residues, dsobj)
        conf_change = score_rmsd(pdb1, pdb2, reslist=reslist)
        rmsd_scores.append((conf_change, conf_change_weight))
        
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
    
    custom_scores = []
    if conf.scoring.custom.score_custom:
        custom_chains = [x.strip() for x in conf.scoring.custom.preds.split(",")]
        custom_weights = parse_weights(conf.scoring.custom.weights, custom_chains)
        i=0
        for result, pred in zip(results, all_chains):
            if pred in custom_chains:
                custom_score = score_custom(result, conf)
                custom_scores.append((custom_score, custom_weights[i]))
                i+=1
    
    other_scores = []
    #TODO add H bond score
    # analyzer = HydrogenBondAnalyzer(pdb_file, distance_cutoff, angle_cutoff)
    # hbonds = analyzer.find_hydrogen_bonds()
    # analyzer.print_hbond_summary(hbonds)
    
    #TODO add hydphob score    
    
    overall_score = 0
    all_scores = plddt_scores + ptm_scores + contacts_scores + rmsd_scores + other_scores + custom_scores
    for s in all_scores:
        overall_score += s[0]*s[1]
    
    all_scores.insert(0, overall_score)
    return all_scores, pdbs, results

def score_custom(result, conf):
    try:
        if conf.scoring.custom:
            scorefile = conf.scoring.custom.custom_script.rsplit("/", 1)
            scorepath = scorefile[0]
            scorefilename = scorefile[1].split(".")[0]

            sys.path.append(scorepath)
            mod = importlib.import_module(scorefilename)
            scorefunc = getattr(mod, conf.scoring.custom.custom_score_func)
            return scorefunc(result, conf)
    except:
        raise ValueError("Invalid score function")