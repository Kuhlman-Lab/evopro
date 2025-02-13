

def parse_input_sequences(d, conf):
    if conf.structure_prediction.structure_prediction_tool == "af2":
        return get_formatted_input_af2(d, conf)
    elif conf.structure_prediction.structure_prediction_tool == "rf2":
        return get_formatted_input_rf2(d, conf)

    
def parse_results(results, conf):
    if conf.structure_prediction.structure_prediction_tool == "af2":
        return parse_results_af2(results, conf)
    
    if conf.structure_prediction.structure_prediction_tool == "rf2":
        return parse_results_rf2(results, conf)
    


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
            
