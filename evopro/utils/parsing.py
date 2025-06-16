import random
import json
import numpy as np
from utils.utils import cif_to_pdb

def parse_input_sequences(d, conf, curr_iter):
    if conf.structure_prediction.structure_prediction_tool == "af2":
        return get_formatted_input_af2(d, conf)
    elif conf.structure_prediction.structure_prediction_tool == "rf2":
        return get_formatted_input_rf2(d, conf)
    elif conf.structure_prediction.structure_prediction_tool == "af3":
        return get_formatted_input_af3(d, conf, curr_iter)
    
def parse_results(results, conf):
    if conf.structure_prediction.structure_prediction_tool == "af2":
        return parse_results_af2(results, conf)
    
    # if conf.structure_prediction.structure_prediction_tool == "rf2":
    #     return parse_results_rf2(results, conf)
    
    if conf.structure_prediction.structure_prediction_tool == "af3":
        return parse_results_af3(results, conf)


def get_formatted_input_af2(d, conf):
    formatted_inputs = []
    for chains in conf.structure_prediction.structure_pred_chains:
        c = list(chains)
        input = [[d.get_chain_sequence(chain) for chain in c]]

        formatted_inputs.append(input)
    return formatted_inputs

def get_formatted_input_rf2(d, conf):
    formatted_inputs = []
    for chains in conf.structure_prediction.structure_pred_chains:
        c = list(chains)
        input = [[d.get_chain_sequence(chain) for chain in c]]

        formatted_inputs.append(input)
    return formatted_inputs

    
def get_formatted_input_af3(d, conf, curr_iter):
    data = []
    custom_msa_dict = {}
    custom_template_dict = {}
    if conf.structure_prediction.custom_msa_chains is not None:
        msa_chains = conf.structure_prediction.custom_msa_chains.split(",")
        msa_paths = conf.structure_prediction.custom_msa_paths.split(",")
        if len(msa_chains) != len(msa_paths):
            msa_paths = [msa_paths[0]] * len(msa_chains)
        for chain, path in zip(msa_chains, msa_paths):
            custom_msa_dict[chain] = path
    if conf.structure_prediction.custom_template_chains is not None:
        template_chains = conf.structure_prediction.custom_template_chains.split(",")
        template_paths = conf.structure_prediction.custom_template_paths.split(",")
        if len(template_chains) != len(template_paths):
            template_paths = [template_paths[0]] * len(template_chains)
        for chain, path in zip(template_chains, template_paths):
            custom_template_dict[chain] = path
    
    for chains in conf.structure_prediction.structure_pred_chains:
        output_dict = {}
        output_dict['name'] = "design"
        if conf.structure_prediction.seed[curr_iter-1] == "random":
            seed = random.randint(0, 1000000000)
            output_dict['modelSeeds'] = [seed]
        else:
            output_dict['modelSeeds'] = [float(conf.structure_prediction.seed[curr_iter-1])]
        c = list(chains)
        for chain in c:
            t = d.get_chain_type(chain)
            seq = d.get_chain_sequence(chain)
            mods = d.chains[chain].modifications
            modifications = []
            for mod in mods:
                if mod["chain"] == chain:
                    modifications.append({"ptmType": mod["type"], "ptmPosition": int(mod["resid"])})
            if "sequences" not in output_dict:  
                output_dict["sequences"] = []
            
            if t == "ligand":
                if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["smiles"]:
                    output_dict["sequences"][t]["id"].append(chain)
                else:
                    seq_dict = {t:{'smiles': seq, "id": [chain]}}
                    if chain in custom_msa_dict:
                        seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                    if chain in custom_template_dict:
                        seq_dict[t]["templates"] = custom_template_dict[chain]

                    output_dict["sequences"].append(seq_dict)
            
            elif t=="rna":
                if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["sequence"]:
                    output_dict["sequences"][t]["id"].append(chain)
                else:
                    seq_dict = {t:{'sequence': parse_rna(seq), "id": [chain]}}
                    if chain in custom_msa_dict:
                        seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                    if chain in custom_template_dict:
                        seq_dict[t]["templates"] = custom_template_dict[chain]

                    output_dict["sequences"].append(seq_dict)
            elif t=="dna":
                if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["sequence"]:
                    output_dict["sequences"][t]["id"].append(chain)
                else:
                    seq_dict = {t:{'sequence': seq.upper(), "id": [chain]}}
                    if chain in custom_msa_dict:
                        seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                    if chain in custom_template_dict:
                        seq_dict[t]["templates"] = custom_template_dict[chain]

                    output_dict["sequences"].append(seq_dict)
            
            else:
                if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["sequence"]:
                    output_dict["sequences"][t]["id"].append(chain)

                else:
                    if len(modifications) > 0:
                        seq_dict = {t:{'sequence': seq, "modifications": modifications, "id": [chain]}}
                    else:
                        seq_dict = {t:{'sequence': seq, "id": [chain]}}
                    if chain in custom_msa_dict:
                        seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                    if chain in custom_template_dict:
                        seq_dict[t]["templates"] = custom_template_dict[chain]

                    output_dict["sequences"].append(seq_dict)
        
        data.append(json.dumps(output_dict))
    
    with open(conf.flags.run_dir + "outputs/af3_data.json", "w") as f:
        f.write(json.dumps(output_dict, indent=4))
    return data

def parse_rna(seq):
    new_seq = ""
    for char in seq:
        if char == "b":
            new_seq += "A"
        elif char == "d":
            new_seq += "C"
        elif char == "h":
            new_seq += "G"
        elif char == "u":
            new_seq += "U"
    
    return new_seq

def parse_results_af2(results, conf):
    from alphafold.common import protein
    parsed_results = []
    pdbs = []
    for result in results:
        parsed_result = {}
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        parsed_result['pdb'] = pdb
        parsed_result['plddt'] = result['plddt']
        parsed_result['pae'] = result['pae_output'][0]
        parsed_result['average_pae'] = result['pae_output'][1]
        parsed_result['ptm'] = result['ptm']
        parsed_result['ranking_confidence'] = result['ranking_confidence']
        parsed_result['structure_module'] = result['structure_module']
        
        if 'iptm' in result:
            parsed_result['iptm'] = result['iptm']

        parsed_results.append(parsed_result)
        
    return parsed_results, pdbs

def parse_results_af3(results, conf):
    
    pdbs = []
    parsed_results = []
    for result in results:
        output_dict = {}

        cif = result['cif_str']
        pdb = cif_to_pdb(cif)
        pdbs.append(pdb)
        output_dict['pdb'] = pdb
        output_dict['token_chain_ids'] = result['token_chain_ids']
        output_dict['token_res_ids'] = result['token_res_ids']
        output_dict['atom_chain_ids'] = result['atom_chain_ids']
        output_dict['atom_plddts'] = result['atom_plddts']
        chains = list(set(output_dict['atom_chain_ids'])) # unique chains
        chains.sort() # sort chains alphabetically
        chain_num_atoms = []
        for i in chains:
            chain_num_atoms.append(output_dict['atom_chain_ids'].count(i))
        
        # get mean pLDDT per chain and overall pLLDT
        overall_pLDDT = 0
        mean_pLDDTs = []

        for j in range(len(chains)):
            sum = 0
            for i in range(len(output_dict['atom_plddts'])):
                # only sum if in matching chain
                if output_dict["atom_chain_ids"][i] == chains[j]:
                    sum = sum + output_dict['atom_plddts'][i]
            # add mean pLDDT for each chain     
            mean_pLDDTs.append(sum/chain_num_atoms[j])
            overall_pLDDT = overall_pLDDT + sum
        
        # calc average
        overall_pLDDT = overall_pLDDT/len(output_dict['atom_plddts'])
        
        # residue_plddt_df = calculate_residue_plddt(pdb, output_dict['atom_plddts'])
        # # protein_plddt = residue_plddt_df[residue_plddt_df['residue_type'] == 'Protein']
        # # hetero_plddt = residue_plddt_df[residue_plddt_df['residue_type'] == 'Hetero']
        # # print("Protein pLDDT", protein_plddt)
        
        # outputs
        output_dict['chains'] = chains
        output_dict['chain_num_atoms'] = chain_num_atoms
        output_dict['mean_pLDDTs'] = mean_pLDDTs
        output_dict['overall_pLDDT'] = overall_pLDDT
        # output_dict['residue_plddt_df'] = residue_plddt_df
        
        if 'iptm' in result:
            output_dict['iptm'] = result['iptm']
        if 'chain_iptm' in result:
            output_dict['chain_iptm'] = result['chain_iptm']
        if 'chain_ptm' in result:
            output_dict['chain_ptm'] = result['chain_ptm']
        if 'chain_pair_iptm' in result:
            output_dict['chain_pair_iptm'] = result['chain_pair_iptm']
        if 'chain_pair_pae_min' in result:
            output_dict['chain_pair_pae_min'] = result['chain_pair_pae_min']

        output_dict['ptm'] = result['ptm']
        output_dict['ranking_confidence'] = result['ranking_score']
        output_dict['fraction_disordered'] = result['fraction_disordered']
        output_dict['has_clash'] = result['has_clash']
        output_dict['pae'] = np.array(result['pae'])
        
        parsed_results.append(output_dict)

    return parsed_results, pdbs
