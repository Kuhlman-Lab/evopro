import sys, os
import string
from Bio.PDB import PDBParser, PDBIO
import time
#from Bio.PDB.Chain import Chain
sys.path.append("/proj/kuhl_lab/evopro")
#from geneticalg_helpers import mutate_by_protein_mpnn
from evopro.genetic_alg.DesignSeq import DesignSeq
from evopro.utils.pdb_parser import get_coordinates_pdb, change_chainid_pdb, append_pdbs
from evopro.genetic_alg.geneticalg_helpers import mutate_by_protein_mpnn

sys.path.append('/proj/kuhl_lab/RFdiffusion/scripts/')
from run_inference_dev import main as run_diffusion_func

chain_names = list(string.ascii_uppercase)

def run_partial_diffusion(config_path, pdb_str, output_folder, contigs=None, provide_seq=None, num_pdbs=10, translate=False):
    start = time.time()
    pdbs = []
    pdb_path = output_folder + "diff.pdb"
    with open(pdb_path, "w") as f:
        f.write(pdb_str)
        
    print("Running partial diffusion on input structure", "random translate:", translate)
    while len(pdbs)<num_pdbs:
        print(contigs, provide_seq)
        pdb, _, _, _ = run_diffusion_func(config_path, 
                                          pdb_path, 
                                          contigs=contigs, 
                                          provide_seq=provide_seq, 
                                          num_designs=num_pdbs, 
                                          design_run=True, 
                                          translate=translate)
        #RENUMBER HERE and write out files for MPNN
        for p in pdb:
            pdbs.append(p)
        
    end = time.time()
    print("Time taken for partial diffusion: {} seconds".format(end-start))
    return pdbs
    
def create_new_seqs_diff_mpnn(startseqs, scored_seqs, num_seqs, run_dir, iter_num, all_seqs = [], af2_preds=["AB", "B"],
                         mpnn_temp="0.1", mpnn_version="s_48_020", mpnn_chains=None, bias_AA_dict=None, bias_by_res_dict=None, 
                         config_path=None, contigs=None, provide_seq=None, diff_random_translate=False):
    
    pool = startseqs.copy()
    
    #create a directory within running dir to run diff and protein mpnn
    output_folder = run_dir + "MPNN_" + str(iter_num) + "/"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    
    rep_obj = startseqs[0]
    if mpnn_chains:
        key_seq = rep_obj.get_sequence_string()
        pdb = None
        for i, pred in zip(range(len(af2_preds)), af2_preds):
            if pred in mpnn_chains:
                if not pdb:
                    pdb = scored_seqs[key_seq]["data"][0]["pdb"][i]
                else:
                    chains, _, _ = get_coordinates_pdb(pdb)
                    new_pdb = scored_seqs[key_seq]["data"][0]["pdb"][i]
                    chains_new, _, _ = get_coordinates_pdb(new_pdb)
                    chains_mod = []
                    chain_ind = 0
                    while len(chains_mod)<len(chains_new):
                        if chain_names[chain_ind] not in chains:
                                chains_mod.append(chain_names[chain_ind])
                        chain_ind+=1

                    for cm, cn in zip(chains_mod, chains_new):
                        new_pdb = change_chainid_pdb(new_pdb, old_chain=cn, new_chain=cm)
                    
                    pdb = append_pdbs(pdb, new_pdb)
                            
    else:
        pdb = scored_seqs[key_seq]["data"][0]["pdb"][0]
    
    if config_path:
        print("Running partial diffusion on best structure before MPNN")
        #check contigs here
        contigs = ['66-66/0 B1-394']
        diff_pdbs = run_partial_diffusion(config_path, 
                                          pdb, 
                                          output_folder, 
                                          contigs=contigs, 
                                          provide_seq=provide_seq, 
                                          num_pdbs=num_seqs-len(pool),
                                          translate=diff_random_translate)
    else:
        raise ValueError("config_path must be specified for partial diffusion!")
    
    pdb_dirs = []
    
    for diff_pdb, j in zip(diff_pdbs, range(len(diff_pdbs))):
        
        pdb_dir = output_folder + "seq_" + str(j) + "/"
        if not os.path.isdir(pdb_dir):
            os.makedirs(pdb_dir)
        with open(output_folder + "seq_" + str(j) + "/seq_" + str(j) + ".pdb", "w") as pdbf:
            pdbf.write(str(diff_pdb))
        pdb_dirs.append(pdb_dir)

    k=0
    inf = 0
    while len(pool) < num_seqs:
        results = mutate_by_protein_mpnn(pdb_dirs[k], rep_obj, mpnn_temp, mpnn_version=mpnn_version, bias_AA_dict=bias_AA_dict, bias_by_res_dict=bias_by_res_dict)
        for result in results:
            inf +=1

            seq = result[-1][-1].strip().split("/")
            newseq_sequence = "".join(seq)
            newseq_sequence_check = ",".join(seq)
            newseqobj = DesignSeq(seq=newseq_sequence, sequence=rep_obj.sequence, mutable=rep_obj.mutable, symmetric=rep_obj.symmetric)

            if newseq_sequence_check not in all_seqs:
                pool.append(newseqobj)
            if inf>=50:
                print("too many mpnn runs without generating a new sequence, using random mutation")
                newseq = rep_obj.mutate()
                if ",".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]]) not in all_seqs:
                    pool.append(newseq)
        k+=1
        if k>=len(startseqs):
            k=0

    return pool


if __name__ == "__main__":
    
    config_path = "/work/users/a/m/amritan/evopro_tests/diff/binder/partial.yaml"
    #pdb_path = "/work/users/a/m/amritan/evopro_tests/diff/binder/design_chain.pdb"
    pdb_path = "/work/users/a/m/amritan/evopro_tests/diff/binder/diff.pdb"
    #pdb_path = "/work/users/a/m/amritan/evopro_tests/diff/binder/ang3_binder.pdb"
    output_folder = "/work/users/a/m/amritan/evopro_tests/diff/binder/"
    #contigs = ['234-234/0 70-70']
    #contigs = ['70-70/0 A1-234']
    contigs = ['70-70/0 B1-234']
    #contigs = ['A1-234/0 70-70']
    #contigs = ['80-80/0 80-80/0 80-80/ 64-64']
    #provide_seq = ['80-112']
    provide_seq=None
    with open(pdb_path, "r") as f:
        pdb_str = f.read()
    pdbs = run_partial_diffusion(config_path, pdb_str, output_folder, contigs=contigs, provide_seq=provide_seq, num_pdbs=3)
    print(pdbs)
    print(len(pdbs))
    with open(output_folder + "temp.pdb", "w") as f:
        f.write(pdb_str)