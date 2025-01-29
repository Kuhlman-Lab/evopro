import numpy as np
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from typing import Tuple, List, Dict, Set
import copy
from dataclasses import dataclass
import logging
from collections import defaultdict

def extract_ca_coordinates(structure: PandasPdb, chain_id: str) -> np.ndarray:
    """
    Extract CA atom coordinates for a specific chain.
    
    Args:
        structure: BioPandas PDB structure
        chain_id: Chain identifier
        
    Returns:
        Array of shape [num_residues, 3] with CA coordinates
    """
    ca_atoms = structure.df['ATOM'][
        (structure.df['ATOM']['chain_id'] == chain_id) &
        (structure.df['ATOM']['atom_name'] == 'CA')
    ]
    
    return ca_atoms[['x_coord', 'y_coord', 'z_coord']].values

def calculate_chain_rmsd(chain1_coords: np.ndarray, chain2_coords: np.ndarray) -> float:
    """
    Calculate RMSD between two chains after centering.
    
    Args:
        chain1_coords: First chain coordinates [num_residues, 3]
        chain2_coords: Second chain coordinates [num_residues, 3]
        
    Returns:
        RMSD value
    """
    if len(chain1_coords) == 0 or len(chain2_coords) == 0:
        return float('inf')
        
    # Take the minimum length
    min_length = min(len(chain1_coords), len(chain2_coords))
    chain1_coords = chain1_coords[:min_length]
    chain2_coords = chain2_coords[:min_length]
    
    # Center coordinates
    chain1_center = chain1_coords.mean(axis=0)
    chain2_center = chain2_coords.mean(axis=0)
    chain1_centered = chain1_coords - chain1_center
    chain2_centered = chain2_coords - chain2_center
    
    # Calculate RMSD
    rmsd = np.sqrt(np.mean(np.sum((chain1_centered - chain2_centered) ** 2, axis=1)))
    
    # Penalize length differences
    length_diff = abs(len(chain1_coords) - len(chain2_coords))
    length_penalty = length_diff / max(len(chain1_coords), len(chain2_coords))
    
    return rmsd * (1 + length_penalty)

def build_chain_distance_matrix(
    pred_structure: PandasPdb,
    true_structure: PandasPdb
) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Build distance matrix between all chains from predicted and true structures.
    
    Args:
        pred_structure: Predicted structure
        true_structure: True structure
        
    Returns:
        Tuple of (distance matrix, pred chain IDs, true chain IDs)
    """
    pred_chains = list(pred_structure.df['ATOM']['chain_id'].unique())
    true_chains = list(true_structure.df['ATOM']['chain_id'].unique())
    
    distances = np.zeros((len(pred_chains), len(true_chains)))
    
    for i, pred_chain in enumerate(pred_chains):
        pred_coords = extract_ca_coordinates(pred_structure, pred_chain)
        for j, true_chain in enumerate(true_chains):
            true_coords = extract_ca_coordinates(true_structure, true_chain)
            distances[i, j] = calculate_chain_rmsd(pred_coords, true_coords)
            
    return distances, pred_chains, true_chains

def get_optimal_chain_mapping(
    pred_structure: PandasPdb,
    true_structure: PandasPdb,
    max_rmsd: float = float('inf')
) -> Dict[str, str]:
    """
    Find optimal mapping between predicted and true chains using Hungarian algorithm.
    
    Args:
        pred_structure: Predicted structure
        true_structure: True structure
        max_rmsd: Maximum RMSD threshold for chain matching
        
    Returns:
        Dictionary mapping predicted chain IDs to true chain IDs
    """
    # Build distance matrix
    distances, pred_chains, true_chains = build_chain_distance_matrix(
        pred_structure, true_structure
    )
    
    # Apply RMSD threshold
    distances[distances > max_rmsd] = float('inf')
    
    # Find optimal assignment
    row_ind, col_ind = linear_sum_assignment(distances)
    
    # Create mapping
    chain_mapping = {}
    for i, j in zip(row_ind, col_ind):
        if distances[i, j] != float('inf'):
            chain_mapping[pred_chains[i]] = true_chains[j]
            
    return chain_mapping

def analyze_chain_alignment(
    pred_structure: PandasPdb,
    true_structure: PandasPdb,
    chain_mapping: Dict[str, str]
) -> Dict[str, Dict[str, float]]:
    """
    Analyze the quality of chain alignment.
    
    Args:
        pred_structure: Predicted structure
        true_structure: True structure
        chain_mapping: Dictionary mapping predicted to true chain IDs
        
    Returns:
        Dictionary of alignment metrics per chain
    """
    metrics = {}
    
    for pred_chain, true_chain in chain_mapping.items():
        pred_coords = extract_ca_coordinates(pred_structure, pred_chain)
        true_coords = extract_ca_coordinates(true_structure, true_chain)
        
        min_length = min(len(pred_coords), len(true_coords))
        rmsd = calculate_chain_rmsd(pred_coords, true_coords)
        
        metrics[pred_chain] = {
            'matched_to': true_chain,
            'rmsd': rmsd,
            'num_residues': len(pred_coords),
            'aligned_residues': min_length,
            'coverage': min_length / max(len(pred_coords), len(true_coords)),
            'sequence_identity': None  # Could be added if sequence information is needed
        }
    
    return metrics

def align_structures(
    pred_structure: PandasPdb,
    true_structure: PandasPdb,
    max_rmsd: float = float('inf')
) -> Tuple[Dict[str, str], Dict[str, Dict[str, float]]]:
    """
    Main function to find optimal chain permutation alignment.
    
    Args:
        pred_structure: Predicted structure
        true_structure: True structure
        max_rmsd: Maximum RMSD threshold for chain matching
        
    Returns:
        Tuple of (chain mapping dictionary, alignment metrics dictionary)
    """
    # Get optimal chain mapping
    chain_mapping = get_optimal_chain_mapping(pred_structure, true_structure, max_rmsd)
    #print(chain_mapping)
    
    # Analyze alignment quality
    metrics = analyze_chain_alignment(pred_structure, true_structure, chain_mapping)
    
    return chain_mapping, metrics

def find_optimal_transformation(
    coords1: np.ndarray,
    coords2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find optimal rotation and translation that minimizes RMSD between two coordinate sets.
    Implementation of Kabsch algorithm.
    
    Args:
        coords1: First set of coordinates [n_points, 3]
        coords2: Second set of coordinates [n_points, 3]
        
    Returns:
        Tuple of (rotation matrix [3,3], translation vector [3])
    """
    # Center the coordinates
    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)
    
    coords1_centered = coords1 - centroid1
    coords2_centered = coords2 - centroid2
    
    # Calculate covariance matrix
    covariance = coords1_centered.T @ coords2_centered
    
    # Singular value decomposition
    U, _, Vt = np.linalg.svd(covariance)
    
    # Handle reflection case
    d = np.linalg.det(Vt.T @ U.T)
    if d < 0:
        Vt[-1] *= -1
    
    # Calculate rotation matrix
    rotation = Vt.T @ U.T
    
    # Calculate translation
    translation = centroid2 - centroid1 @ rotation
    
    return rotation, translation

def transform_structure(
    structure: PandasPdb,
    rotation: np.ndarray,
    translation: np.ndarray
) -> PandasPdb:
    """
    Apply rotation and translation to a structure.
    
    Args:
        structure: PDB structure
        rotation: 3x3 rotation matrix
        translation: Translation vector
        
    Returns:
        Transformed structure
    """
    new_structure = copy.copy(structure)
    coords = new_structure.df['ATOM'][['x_coord', 'y_coord', 'z_coord']].values
    
    # Apply transformation
    new_coords = coords @ rotation + translation
    
    new_structure.df['ATOM'][['x_coord', 'y_coord', 'z_coord']] = new_coords
    return new_structure

def align_and_save_structure(
    pred_structure: PandasPdb,
    true_structure: PandasPdb,
    chain_mapping: Dict[str, str],
    output_path: str = None,
    max_rmsd: float = float('inf')
) -> Tuple[PandasPdb, float]:
    """
    Find best superposition based on chain mapping and save aligned structure.
    
    Args:
        pred_structure: Predicted structure
        true_structure: True structure
        chain_mapping: Dictionary mapping predicted to true chain IDs
        output_path: Path to save aligned structure
        max_rmsd: Maximum RMSD threshold for including chains in alignment
        
    Returns:
        Tuple of (transformed structure, final RMSD)
    """
    # Collect all CA atoms from matched chains
    pred_cas = []
    true_cas = []
    
    for pred_chain, true_chain in chain_mapping.items():
        pred_coords = extract_ca_coordinates(pred_structure, pred_chain)
        true_coords = extract_ca_coordinates(true_structure, true_chain)
        
        # Take minimum length for each chain pair
        min_length = min(len(pred_coords), len(true_coords))
        
        if min_length > 0:
            # Calculate initial RMSD for this chain pair
            initial_rmsd = calculate_chain_rmsd(pred_coords[:min_length], true_coords[:min_length])
            
            # Only include chains with RMSD below threshold
            if initial_rmsd <= max_rmsd:
                pred_cas.append(pred_coords[:min_length])
                true_cas.append(true_coords[:min_length])
    
    if not pred_cas:
        raise ValueError("No chains found within RMSD threshold for alignment")
    
    # Concatenate coordinates from all matched chains
    pred_cas = np.vstack(pred_cas)
    true_cas = np.vstack(true_cas)
    
    # Find optimal transformation
    rotation, translation = find_optimal_transformation(pred_cas, true_cas)
    #print(rotation, translation)
    # Apply transformation to entire structure
    aligned_structure = transform_structure(pred_structure, rotation, translation)
    #print(aligned_structure)
    # Calculate final RMSD
    transformed_coords = pred_cas @ rotation + translation
    final_rmsd = np.sqrt(np.mean(np.sum((transformed_coords - true_cas) ** 2, axis=1)))
    
    # Save aligned structure
    if output_path is not None:
        aligned_structure.to_pdb(output_path)
    
    return aligned_structure, final_rmsd

def align_structures_and_save(
    pred_structure: PandasPdb,
    true_structure: PandasPdb,
    output_path: str = None,
    max_rmsd: float = float('inf')
) -> Tuple[Dict[str, str], Dict[str, Dict[str, float]], float]:
    """
    Complete pipeline: find optimal chain mapping, align structure, and save PDB.
    
    Args:
        pred_structure: Predicted structure
        true_structure: True structure
        output_path: Path to save aligned structure
        max_rmsd: Maximum RMSD threshold for chain matching
        
    Returns:
        Tuple of (chain mapping, alignment metrics, final RMSD)
    """
    # Get optimal chain mapping
    chain_mapping = get_optimal_chain_mapping(pred_structure, true_structure, max_rmsd)
    
    # Align structure and save PDB
    aligned_structure, final_rmsd = align_and_save_structure(
        pred_structure, true_structure, chain_mapping, output_path=output_path, max_rmsd=max_rmsd
    )
    
    # Get alignment metrics for transformed structure
    metrics = analyze_chain_alignment(aligned_structure, true_structure, chain_mapping)
    
    return chain_mapping, metrics, final_rmsd

if __name__ == "__main__":
    
    #THIS DOESNT WORK CORRECTLY! SEE multichain_permutation_align_test.py
    path = "/work/users/a/m/amritan/evopro_tests/rmsd/multichain_perm/test6/"
    pdb1_path = path + "seq_0_iter_1_model_1_chainABCDE_round3_badRMSD.pdb"
    pdb2_path = path + "rfdiff_model25step_0_renumbered_startmodel.pdb"
    output_path = path + "aligned_badRMSD.pdb"
    # output_path = None
    
    pred_struct = PandasPdb().read_pdb(pdb1_path)
    true_struct = PandasPdb().read_pdb(pdb2_path)
    
    try:
        # Find optimal chain mapping, align structure, and save PDB
        chain_mapping, metrics, final_rmsd = align_structures_and_save(
            pred_struct,
            true_struct,
            output_path = output_path,
            max_rmsd=50.0
        )
        
        print(f"RMSD: {final_rmsd:.2f} Å")
        # print("\nChain Mapping Results:")
        # for pred_chain, true_chain in chain_mapping.items():
        #     print(f"\nPredicted chain {pred_chain} matched to true chain {true_chain}:")
        #     m = metrics[pred_chain]
        #     print(f"  RMSD: {m['rmsd']:.2f} Å")
        #     print(f"  Residues: {m['num_residues']}")
        #     print(f"  Aligned residues: {m['aligned_residues']}")
        #     print(f"  Coverage: {m['coverage']*100:.1f}%")
        
        # # Print unmatched chains
        # pred_unmatched = set(pred_struct.df['ATOM']['chain_id'].unique()) - set(chain_mapping.keys())
        # true_unmatched = set(true_struct.df['ATOM']['chain_id'].unique()) - set(chain_mapping.values())
        
        # if pred_unmatched:
        #     print("\nUnmatched predicted chains:", ", ".join(pred_unmatched))
        # if true_unmatched:
        #     print("Unmatched true chains:", ", ".join(true_unmatched))
            
    except Exception as e:
        print(f"Error during alignment: {e}")
        
    path = "/work/users/a/m/amritan/evopro_tests/rmsd/multichain_perm/test6/"
    pdb1_path = path + "seq_0_final_model_1_chainABCDE_round2_goodRMSD.pdb"
    pdb2_path = path + "rfdiff_model25step_0_renumbered_startmodel.pdb"
    output_path = path + "aligned_goodRMSD.pdb"
    # output_path = None
    
    pred_struct = PandasPdb().read_pdb(pdb1_path)
    true_struct = PandasPdb().read_pdb(pdb2_path)
    
    try:
        # Find optimal chain mapping, align structure, and save PDB
        chain_mapping, metrics, final_rmsd = align_structures_and_save(
            pred_struct,
            true_struct,
            output_path = output_path,
            max_rmsd=50.0
        )
        
        print(f"RMSD: {final_rmsd:.2f} Å")
        # print("\nChain Mapping Results:")
        # for pred_chain, true_chain in chain_mapping.items():
        #     print(f"\nPredicted chain {pred_chain} matched to true chain {true_chain}:")
        #     m = metrics[pred_chain]
        #     print(f"  RMSD: {m['rmsd']:.2f} Å")
        #     print(f"  Residues: {m['num_residues']}")
        #     print(f"  Aligned residues: {m['aligned_residues']}")
        #     print(f"  Coverage: {m['coverage']*100:.1f}%")
        
        # # Print unmatched chains
        # pred_unmatched = set(pred_struct.df['ATOM']['chain_id'].unique()) - set(chain_mapping.keys())
        # true_unmatched = set(true_struct.df['ATOM']['chain_id'].unique()) - set(chain_mapping.values())
        
        # if pred_unmatched:
        #     print("\nUnmatched predicted chains:", ", ".join(pred_unmatched))
        # if true_unmatched:
        #     print("Unmatched true chains:", ", ".join(true_unmatched))
            
    except Exception as e:
        print(f"Error during alignment: {e}")