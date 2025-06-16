import gzip
import io
import copy
from typing import Any, List, Optional, Union, Dict
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial.distance import cdist

import utils.hbond_constants as rc
from utils.utils import build_sc_from_chi, impute_CB, get_worker_rng
from objects.graph import Graph, GraphLike

# Complete sequence of chain IDs supported by the PDB format.
PDB_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
PDB_MAX_CHAINS = len(PDB_CHAIN_IDS)  # := 62.


class Protein(GraphLike):
    """Protein structure representation."""

    # Cartesian coordinates of atoms in angstroms. The atom types correspond to
    # residue_constants.atom_types, i.e. the first three are N, CA, C.
    atom27_xyz: np.ndarray  # [num_res, 27, 3]

    # Amino-acid type for each residue represented as an integer between 0 and
    # 20, where 20 is 'X'.
    aatype: np.ndarray  # [num_res]

    # Binary float mask to indicate presence of a particular atom. 1.0 if an atom
    # is present and 0.0 if not. This should be used for loss masking.
    atom27_mask: np.ndarray  # [num_res, 27]

    # Residue index as used in PDB. It is not necessarily continuous or 0-indexed.
    residue_index: np.ndarray  # [num_res]

    # 0-indexed number corresponding to the chain in the protein that this residue
    # belongs to.
    chain_index: np.ndarray  # [num_res]

    # B-factors, or temperature factors, of each residue (in sq. angstroms units),
    # representing the displacement of the residue from its ground truth mean
    # value.
    b_factors: np.ndarray  # [num_res, 27]

    # Boolean representing whether a residue is "designable" or not. Designable
    # doesn't necessary refer to the amino acid sequence; A residue is designable
    # if it can be modified (e.g. packing the sidechain).
    designable_res: np.ndarray  # [num_res]

    # Boolean representing which residue(s) are "neighbors" of each other. Used to narrow down search space.
    neighbor_mask: np.ndarray  # [num_res, num_res]

    # Dictionary containing properties of HETATM associated with the protein. Properties include
    # atom_name, element, res_name, residue_index, chain_index, and atom_xyz.
    hetatm_dict: Dict[str, np.ndarray] = {}

    def __init__(
        self,
        atom27_xyz: np.ndarray,
        aatype: np.ndarray,
        atom27_mask: np.ndarray,
        residue_index: np.ndarray,
        chain_index: np.ndarray,
        b_factors: np.ndarray,
        designable_res: Optional[np.ndarray] = None,
        neighbor_mask: Optional[np.ndarray] = None,
        hetatm_dict: Dict[str, np.ndarray] = {},
        graph: Optional[Graph] = None,
    ) -> None:
        # Protein properties
        self.atom27_xyz = atom27_xyz
        self.aatype = aatype
        self.atom27_mask = atom27_mask
        self.residue_index = residue_index
        self.chain_index = chain_index
        self.b_factors = b_factors
        self.hetatm_dict = hetatm_dict

        if designable_res is None:
            self.designable_res = np.ones_like(self.aatype).astype(bool)
        else:
            self.designable_res = designable_res

        if neighbor_mask is None:
            self.neighbor_mask = np.ones((self.aatype.size, self.aatype.size)).astype(
                bool
            )
        else:
            self.neighbor_mask = neighbor_mask

        self._validate_shapes()
        if len(np.unique(self.chain_index)) > PDB_MAX_CHAINS:
            raise ValueError(
                f"Cannot build an instance with more than {PDB_MAX_CHAINS} chains "
                "because these cannot be written to PDB format."
            )

        # Initialize the graph
        if graph is None:
            self.reset_graph()
        else:
            self.graph = graph

    def _validate_shapes(self) -> None:
        self.n_res = self.atom27_xyz.shape[0]

        # If atom14 arrays were provided, pad them to atom27.
        if (
            self.atom27_xyz.shape[1] == 14
            and self.atom27_mask.shape[1] == 14
            and self.b_factors.shape[1] == 14
        ):
            self.atom27_xyz = np.concatenate(
                [self.atom27_xyz, np.zeros((self.n_res, 13, 3))], axis=1
            )
            self.atom27_mask = np.concatenate(
                [self.atom27_mask, np.zeros((self.n_res, 13))], axis=1
            )
            self.b_factors = np.concatenate(
                [self.b_factors, np.zeros((self.n_res, 13))], axis=1
            )

        # Validate the shapes of the arrays.
        assert self.atom27_xyz.shape == (self.n_res, 27, 3)
        assert self.aatype.shape == (self.n_res,)
        assert self.atom27_mask.shape == (self.n_res, 27)
        assert self.residue_index.shape == (self.n_res,)
        assert self.chain_index.shape == (self.n_res,)
        assert self.b_factors.shape == (self.n_res, 27)
        assert self.designable_res.shape == (self.n_res,)
        assert self.neighbor_mask.shape == (self.n_res, self.n_res)

        # Validate hetatm_dict if present.
        if self.hetatm_dict != {}:
            assert "atom_name" in self.hetatm_dict
            assert "element" in self.hetatm_dict
            assert "res_name" in self.hetatm_dict
            assert "residue_index" in self.hetatm_dict
            assert "chain_index" in self.hetatm_dict
            assert "atom_xyz" in self.hetatm_dict
            assert (
                self.hetatm_dict["atom_name"].shape[0]
                == self.hetatm_dict["element"].shape[0]
            )
            assert (
                self.hetatm_dict["atom_name"].shape[0]
                == self.hetatm_dict["res_name"].shape[0]
            )
            assert (
                self.hetatm_dict["atom_name"].shape[0]
                == self.hetatm_dict["residue_index"].shape[0]
            )
            assert (
                self.hetatm_dict["atom_name"].shape[0]
                == self.hetatm_dict["chain_index"].shape[0]
            )
            assert (
                self.hetatm_dict["atom_name"].shape[0]
                == self.hetatm_dict["atom_xyz"].shape[0]
            )

    @classmethod
    def from_pdb_string(
        cls,
        pdb_str: str,
        model_idx: int = 0,
        chain_id: Optional[Union[str, List[str]]] = None,
        discard_water: bool = True,
        discard_Hs: bool = True,
        mse_to_met: bool = False,
        keep_hetatm: bool = False,
    ):
        """Takes a PDB string and constructs an object of class cls.

        WARNING: All non-standard residue types will be converted into UNK. All
            non-standard atoms will be ignored.

        Args:
            cls: Class of resulting object.
            pdb_str (str): The contents of the PDB file
            model_idx (int, optional): The specific model in the PDB file that will be
                parsed. This is 0-indexed. Defaults to 0.
            chain_id (Union[str, List[str]], optional): If chain_id is specified
                (e.g. ['A']), then only those chains are parsed. Otherwise all chains
                are parsed. Defaults to None.
            discard_water (bool, optional): Boolean specifying whether to ignore water molecules.
                Defaults to True.
            discard_Hs (bool, optional): Boolean specifying whether to ignore hydrogen atoms.
                If hydrogens are discarded, then the residues will be reduced to 14 atoms per
                residue. Otherwise, the residues will be have up to 27 atoms per residue.
                Defaults to True.
            mse_to_met (bool, optional): Boolean specifying whether to convert MSE residues to
                MET residues. Defaults to False.
            keep_hetatm (bool, optional): Boolean specifying whether to keep HETATM records. If
                True, then they will be added to the hetatm_dict. Defaults to False.

        Returns:
            A new cls object parsed from the PDB contents.
        """
        pdb_fh = io.StringIO(pdb_str)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("none", pdb_fh)
        models = list(structure.get_models())
        if model_idx is not None and model_idx > len(models) - 1:
            raise ValueError(
                f"Requested model index is out of range. Found {len(models)} models."
            )
        elif model_idx is not None:
            model = models[model_idx]
        else:
            model = models[0]
        if isinstance(chain_id, str):
            chain_id = [chain_id]

        if discard_Hs:
            n_atoms = 14
            restype_name_to_atom_name = rc.restype_name_to_atom14_names
        else:
            n_atoms = 27
            restype_name_to_atom_name = rc.restype_name_to_atom27_names

        atom_positions = []
        aatype = []
        atom_mask = []
        residue_index = []
        chain_ids = []
        b_factors = []
        hetatm_dict = {
            "atom_name": [],
            "element": [],
            "res_name": [],
            "residue_index": [],
            "chain_index": [],
            "atom_xyz": [],
        }
        insertion_code_offset = 0
        hetatm_insertion_code_offset = 0
        for chain in sorted(model, key=lambda x: x.id):
            if chain_id is not None and chain.id not in chain_id:
                continue
            for res in sorted(chain, key=lambda x: x.id[1]):
                # Discard water residues.
                if discard_water:
                    if res.resname == "HOH":
                        continue

                # Convert MSE residues to MET by changing only necessary fields.
                if mse_to_met:
                    if res.resname == "MSE":
                        res.id = (" ", *res.id[1:])
                        res.resname = "MET"
                        for atom in res:
                            if atom.name == "SE":
                                atom.name = "SD"

                # Parse for hetatms.
                is_het_res = res.id[0] != " "
                if keep_hetatm and is_het_res:
                    if res.id[2] != " ":
                        hetatm_insertion_code_offset += 1

                    atom_name = []
                    atom_elem = []
                    atom_xyz = []
                    for atom in res:
                        atom_name.append(atom.name)
                        atom_elem.append(atom.element)
                        atom_xyz.append(atom.coord)
                    res_idx = [res.id[1] + hetatm_insertion_code_offset] * len(
                        atom_elem
                    )
                    chain_idx = [chain.id] * len(atom_elem)
                    res_name = [res.resname] * len(atom_elem)

                    hetatm_dict["atom_name"].extend(atom_name)
                    hetatm_dict["element"].extend(atom_elem)
                    hetatm_dict["res_name"].extend(res_name)
                    hetatm_dict["residue_index"].extend(res_idx)
                    hetatm_dict["chain_index"].extend(chain_idx)
                    hetatm_dict["atom_xyz"].extend(atom_xyz)
                    continue
                elif not keep_hetatm and is_het_res:
                    continue

                # Increment residue index offset if insertion code is detected.
                if res.id[2] != " ":
                    insertion_code_offset += 1

                res_shortname = rc.restype_3to1.get(res.resname, "X")
                if res_shortname == "X":
                    # Skip non-standard residues that remain.
                    continue
                restype_idx = rc.restype_order.get(res_shortname, rc.restype_num)
                pos = np.zeros((n_atoms, 3))
                mask = np.zeros((n_atoms,))
                res_b_factors = np.zeros((n_atoms,))
                for atom in res:
                    if atom.name not in restype_name_to_atom_name[res.resname]:
                        continue
                    pos[restype_name_to_atom_name[res.resname].index(atom.name)] = (
                        atom.coord
                    )
                    mask[restype_name_to_atom_name[res.resname].index(atom.name)] = 1.0
                    res_b_factors[
                        restype_name_to_atom_name[res.resname].index(atom.name)
                    ] = atom.bfactor
                if np.sum(mask) < 0.5:
                    # If no known atom positions are reported for the residue then skip it.
                    continue

                # Update protein-level lists
                aatype.append(restype_idx)
                atom_positions.append(pos)
                atom_mask.append(mask)
                residue_index.append(res.id[1] + insertion_code_offset)
                chain_ids.append(chain.id)
                b_factors.append(res_b_factors)

        # Chain IDs are usually characters so map these to ints using PDB_CHAIN_IDS.
        unique_chain_ids = np.unique(chain_ids + hetatm_dict["chain_index"])
        chain_id_mapping = {cid: PDB_CHAIN_IDS.index(cid) for cid in unique_chain_ids}
        chain_index = np.array([chain_id_mapping[cid] for cid in chain_ids])

        # Convert hetatm_dict components to arrays.
        if len(hetatm_dict["atom_name"]) > 0:
            hetatm_dict["chain_index"] = [
                chain_id_mapping[cid] for cid in hetatm_dict["chain_index"]
            ]
            hetatm_dict = {k: np.stack(v) for k, v in hetatm_dict.items()}
        else:
            hetatm_dict = {}

        return cls(
            atom27_xyz=np.array(atom_positions),
            atom27_mask=np.array(atom_mask),
            aatype=np.array(aatype),
            residue_index=np.array(residue_index),
            chain_index=chain_index,
            b_factors=np.array(b_factors),
            hetatm_dict=hetatm_dict if keep_hetatm else {},
        )

    @classmethod
    def from_pdb_file(cls, pdb_file: Union[str, io.StringIO], **kwargs):
        if isinstance(pdb_file, str):
            # Obtain PDB string from PDB file.
            if pdb_file[-3:] == "pdb":
                with open(pdb_file, "r") as f:
                    pdb_str = f.read()
            elif pdb_file[-6:] == "pdb.gz":
                with gzip.open(pdb_file, "rb") as f:
                    pdb_str = f.read().decode()
            else:
                raise ValueError("Unrecognized file type.")
        else:
            # Get string from the StringIO
            pdb_str = pdb_file.read()

        # Parse the string and get Protein.
        return cls.from_pdb_string(pdb_str, **kwargs)

    @staticmethod
    def _chain_end(
        atom_index: str, end_resname: str, chain_name: str, residue_index: int
    ) -> str:
        chain_end = "TER"
        return (
            f"{chain_end:<6}{atom_index:>5}      {end_resname:>3} "
            f"{chain_name:>1}{residue_index:>4}"
        )

    def to_pdb(
        self, skip_unk: bool = False, unk_to_gly: bool = False, no_hetatm: bool = True
    ) -> str:
        """Converts `Protein` instance to a PDB string.

        Arguments:
            skip_unk: Skip unknown residues and only output those with defined restypes. Default is False.
            unk_to_gly: Convert unknown residues to glycines. Default is False.
            no_hetatm (bool, optional): Disables outputting HETATM records. Defaults to True.

        Returns:
            PDB string.
        """

        if skip_unk and unk_to_gly:
            raise ValueError(
                "Cannot specify both skip_unk and unk_to_gly at the same time."
            )

        restypes = rc.restypes + ["X"]

        def res_1to3(r):
            return rc.restype_1to3.get(restypes[r], "UNK")

        pdb_lines = []

        atom_mask = self.atom27_mask
        aatype = copy.deepcopy(self.aatype)  # Copy to avoid changing underlying data
        atom_positions = self.atom27_xyz
        residue_index = self.residue_index.astype(np.int32)
        chain_index = self.chain_index.astype(np.int32)
        b_factors = self.b_factors

        if np.any(aatype > rc.restype_num):
            raise ValueError("Invalid aatypes.")

        if unk_to_gly:
            unk_mask = aatype == rc.restype_num
            aatype[unk_mask] = rc.restype_order["G"]

        # Construct a mapping from chain integer indices to chain ID strings.
        chain_ids = {}
        if not no_hetatm and self.hetatm_dict != {}:
            all_chain_index = np.concatenate(
                [chain_index, self.hetatm_dict["chain_index"]]
            )
        else:
            all_chain_index = chain_index
        for i in np.unique(all_chain_index):  # np.unique gives sorted output.
            if i >= PDB_MAX_CHAINS:
                raise ValueError(
                    f"The PDB format supports at most {PDB_MAX_CHAINS} chains."
                )
            chain_ids[i] = PDB_CHAIN_IDS[i]

        pdb_lines.append("MODEL     1")
        atom_index = 1
        last_chain_index = chain_index[0]
        # Add all atom sites.
        for i in range(aatype.shape[0]):
            if skip_unk and (aatype[i] == rc.restype_num):
                continue
            # Close the previous chain if in a multichain PDB.
            if last_chain_index != chain_index[i]:
                pdb_lines.append(
                    self._chain_end(
                        atom_index,
                        res_1to3(aatype[i - 1]),
                        chain_ids[chain_index[i - 1]],
                        residue_index[i - 1],
                    )
                )
                last_chain_index = chain_index[i]
                atom_index += 1  # Atom index increases at the TER symbol.

            res_name_3 = res_1to3(aatype[i])
            atom_types = rc.restype_name_to_atom27_names[res_1to3(aatype[i])]
            for atom_name, pos, mask, b_factor in zip(
                atom_types, atom_positions[i], atom_mask[i], b_factors[i]
            ):
                if mask < 0.5:
                    continue

                record_type = "ATOM"
                name = atom_name if len(atom_name) == 4 else f" {atom_name}"
                alt_loc = ""
                insertion_code = ""
                occupancy = 1.00
                element = atom_name[0]  # Protein supports only C, N, O, S, this works.
                charge = ""

                # PDB is a columnar format, every space matters here!
                atom_line = (
                    f"{record_type:<6}{atom_index:>5} {name:<4}{alt_loc:>1}"
                    f"{res_name_3:>3} {chain_ids[chain_index[i]]:>1}"
                    f"{residue_index[i]:>4}{insertion_code:>1}   "
                    f"{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}"
                    f"{occupancy:>6.2f}{b_factor:>6.2f}          "
                    f"{element:>2}{charge:>2}"
                )
                pdb_lines.append(atom_line)
                atom_index += 1

        # Close the final chain.
        pdb_lines.append(
            self._chain_end(
                atom_index,
                res_1to3(aatype[-1]),
                chain_ids[chain_index[-1]],
                residue_index[-1],
            )
        )

        # Potentially add HETATM records.
        if not no_hetatm:
            if self.hetatm_dict != {}:
                for i in range(self.hetatm_dict["element"].shape[0]):
                    record_type = "HETATM"
                    atom_name = self.hetatm_dict["atom_name"][i]
                    name = atom_name if len(atom_name) == 4 else f" {atom_name}"
                    alt_loc = ""
                    insertion_code = ""
                    occupancy = 1.00
                    element = self.hetatm_dict["element"][i]
                    charge = ""
                    pos = self.hetatm_dict["atom_xyz"][i]

                    # PDB is a columnar format, every space matters here!
                    atom_line = (
                        f"{record_type:<6}{atom_index:>5} {name:<4}{alt_loc:>1}"
                        f"{self.hetatm_dict['res_name'][i]:>3} {chain_ids[self.hetatm_dict['chain_index'][i]]:>1}"
                        f"{self.hetatm_dict['residue_index'][i]:>4}{insertion_code:>1}   "
                        f"{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}"
                        f"{occupancy:>6.2f}{b_factor:>6.2f}          "
                        f"{element:>2}{charge:>2}"
                    )
                    pdb_lines.append(atom_line)
                    atom_index += 1

        # End the model and file.
        pdb_lines.append("ENDMDL")
        pdb_lines.append("END")

        # Pad all lines to 80 characters.
        pdb_lines = [line.ljust(80) for line in pdb_lines]
        return "\n".join(pdb_lines) + "\n"  # Add terminating newline.

    def clear_sequence(self, mask: Optional[np.ndarray] = None) -> None:
        """Clears the sequence information from the protein.

        This sets all residues in the mask to an X and strips sidechain atoms (and Hs).
        """
        if mask is None:
            mask = np.ones_like(self.aatype)
        self.aatype = (
            mask * rc.restype_num * np.ones_like(self.aatype) + (1 - mask) * self.aatype
        )
        self.atom27_xyz[:, 4:] = (1 - mask)[..., None, None] * self.atom27_xyz[:, 4:]
        self.atom27_mask[:, 4:] = (1 - mask)[..., None] * self.atom27_mask[:, 4:]

    def clear_sidechains(self, mask: Optional[np.ndarray] = None) -> None:
        """Clears the sidechain coordinates (and Hs) from all residues in the mask."""
        if mask is None:
            mask = np.ones_like(self.aatype)
        self.atom27_xyz[:, 4:] = (1 - mask)[..., None, None] * self.atom27_xyz[:, 4:]
        self.atom27_mask[:, 4:] = (1 - mask)[..., None] * self.atom27_mask[:, 4:]

    def apply_graph_to_protein(self, seq=True, sc=True) -> None:
        """Applies chi angles and aatypes from graph to protein.
        Needed for motif design data loading."""
        nodes = self.graph.nodes
        n_res = self.aatype.size

        # grab chi angles
        chi_angles = np.zeros((n_res, 4))
        chi_angle_mask = np.zeros((n_res, 4), dtype=bool)
        for i in nodes:
            chi_angles[i] = nodes[i]["chi_angles"]
            chi_angle_mask[i] = nodes[i]["chi_angle_mask"]
            if seq:
                self.aatype[i] = nodes[i]["aatype"]

        if sc:
            atom27_xyz, atom27_mask = build_sc_from_chi(
                self.atom27_xyz[:, :4], self.aatype, chi_angles, chi_angle_mask
            )
            nodes = np.array(nodes)
            self.atom27_xyz[nodes, :14] = atom27_xyz[nodes]
            self.atom27_mask[nodes, :14] = atom27_mask[nodes]

    def unk_to_gly(self) -> None:
        """Convert any 'X' residues to 'G' for compatibility with PDB writer."""
        x_mask = self.aatype == rc.restype_num
        self.aatype[x_mask] = rc.restype_order["G"]

    def set_designable_res(self, residues: List[int]) -> None:
        """Sets designable res based on input list or array of res indices."""
        self.designable_res = np.zeros_like(self.designable_res).astype(bool)
        self.designable_res[residues] = True

    def set_neighbor_mask(self, dist: float = 5.0) -> None:
        """Sets nearest neighbor mask based on Cb distance (Angstrom) cutoff."""
        dmat = self._get_Cb_dmat()
        self.neighbor_mask = dmat < dist

    def __hash__(self) -> int:
        """
        Need to define a hash function to compare different GraphLike objects.
        Generates a unique integer based on the specific data in this instance of GraphLike.
        Only immutable objects are hashable.

        For some Protein objects, the PDB fields contain all the unique info needed.
        Therefore, we can just use the scaffold string representation.

        If you need to distinguish between Proteins with the same PDB but different graphs,
        you will need to hash and sort the graph node attributes as a tuple and add these.

        # graph_data = dict(self.hbnet_graph.nodes)
        # resids, aatypes, chi_angles = [], [], []
        # for key, value in graph_data.items():
        #     resids.append(key)
        #     aatypes.append(value['aatype'])
        #     chi_angles.append(tuple(value['chi_angles']))

        # resids, aatypes, chi_angles = (zip(*sorted(zip(resids, aatypes, chi_angles))))
        # graph_hash = hash((resids, aatypes, chi_angles))
        """
        return hash(self.to_pdb(unk_to_gly=True))

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, self.__class__) and (
            self.__hash__() == other.__hash__()
        )

    def _get_Cb_dmat(
        self,
    ) -> np.ndarray:
        res_xyz_all = self.atom27_xyz
        res_xyz_Cb = impute_CB(
            res_xyz_all[:, 0, :], res_xyz_all[:, 1, :], res_xyz_all[:, 2, :]
        )
        return cdist(res_xyz_Cb, res_xyz_Cb)

    def crop(
        self,
        topk: int = 64,
        pad: bool = False,
        mask: Optional[np.ndarray] = None,
    ):
        """Randomly selects a residue and its topk nearest neighbors in 3D space,
        discarding all remaining residues.

        Args:
            topk (int, optional): Number of nearest neighbors to keep after
                cropping, including the seed residue. Defaults to 64.
            pad (bool, optional): Whether to pad proteins too small to fit the crop
                size. Defaults to False.
            mask (np.ndarray, optional): Mask the specifies which residues can be
                selected as seeds for the crop. Defaults to None.

        Returns:
            self.__class__: Cropped Protein object.
        """

        if self.n_res <= topk:
            # No-padding means just return as-is
            if not pad:
                return self
            # Padding means we need to add empty residues
            else:
                return self.pad(topk)
        else:
            protein = self

        if mask is None:
            # If no mask, pull from all residues with complete backbones
            mask = np.prod(protein.atom27_mask[:, :4], axis=-1).astype(bool)
        else:
            assert mask.shape == protein.residue_index.shape

        # Get the seed residue
        rng = get_worker_rng()
        raw_idx = np.arange(protein.residue_index.size)
        seed_res = rng.choice(raw_idx[mask], size=1)

        # Find KNN from distance matrix
        dmat = np.squeeze(protein._get_Cb_dmat()[seed_res])
        # Include seed residue as one of the topk
        knn = np.argpartition(dmat, topk)[:topk]
        # Sort to avoid scrambling residue order
        knn = np.sort(knn)

        # Making new protein will automatically reset the graph and validate shapes
        return self.__class__(
            atom27_xyz=np.take_along_axis(
                protein.atom27_xyz, indices=knn[:, None, None], axis=0
            ),
            aatype=np.take_along_axis(protein.aatype, indices=knn, axis=0),
            atom27_mask=np.take_along_axis(
                protein.atom27_mask, indices=knn[:, None], axis=0
            ),
            residue_index=np.take_along_axis(
                protein.residue_index, indices=knn, axis=0
            ),
            chain_index=np.take_along_axis(protein.chain_index, indices=knn, axis=0),
            b_factors=np.take_along_axis(
                protein.b_factors, indices=knn[:, None], axis=0
            ),
            hetatm_dict=self.hetatm_dict,
        )

    def mask(self, res_idx: np.ndarray):
        """
        Args:
            res_idx (np.ndarray): Array of residue indices to keep after masking.

        Returns:
            self.__class__: Masked Protein object.
        """
        assert res_idx.size <= self.n_res
        return self.__class__(
            atom27_xyz=np.take_along_axis(
                self.atom27_xyz, indices=res_idx[:, None, None], axis=0
            ),
            atom27_mask=np.take_along_axis(
                self.atom27_mask, indices=res_idx[:, None], axis=0
            ),
            aatype=np.take_along_axis(self.aatype, indices=res_idx, axis=0),
            residue_index=np.take_along_axis(
                self.residue_index, indices=res_idx, axis=0
            ),
            chain_index=np.take_along_axis(self.chain_index, indices=res_idx, axis=0),
            b_factors=np.take_along_axis(
                np.zeros_like(self.atom27_mask), indices=res_idx[:, None], axis=0
            ),
        )

    def crop_contiguous(
        self,
        size: int = 64,
        pad: bool = False,
        mask: Optional[np.ndarray] = None,
    ):
        """Randomly selects a residue and the next size neighbors in seq dimension, discarding the others.

        Args:
            size (int, optional): Total size of crop. Defaults to 64.
            pad (bool, optional): Whether to pad proteins too small to fit the crop
                size. Defaults to False.
            mask (np.ndarray, optional): Mask the specifies which residues can be
                selected as seeds for the crop. Defaults to None.
        Returns:
            self.__class__: Cropped Protein object.
        """

        if self.n_res <= size:
            # No-padding means just return as-is
            if not pad:
                return self
            # Padding means we need to add empty residues
            else:
                return self.pad(size)
        else:
            protein = self

            # Get the seed residue
            rng = get_worker_rng()
            raw_idx = np.arange(protein.residue_index.size)
            seed_res = rng.choice(raw_idx[mask], size=1)[0]

            # Check if crop will overflow length of the protein
            overflow = (seed_res + size) - protein.n_res
            if overflow > 0:
                # If so, pad the protein
                protein = protein.pad(protein.n_res + overflow)

        # Making new protein will automatically reset the graph and validate shapes
        return self.__class__(
            atom27_xyz=protein.atom27_xyz[seed_res : seed_res + size],
            aatype=protein.aatype[seed_res : seed_res + size],
            atom27_mask=protein.atom27_mask[seed_res : seed_res + size],
            residue_index=protein.residue_index[seed_res : seed_res + size],
            chain_index=protein.chain_index[seed_res : seed_res + size],
            b_factors=protein.b_factors[seed_res : seed_res + size],
            hetatm_dict=self.hetatm_dict,
        )

    def pad(self, n: int = 128):
        """
        Pad Protein to an arbitrary length by appending ghost residues to the end.

        """
        n_atoms = self.atom27_xyz.shape[1]
        n_add = n - self.n_res
        if n_add < 0:
            raise ValueError(
                "Padded size {n} is smaller than actual size {self.n_res}."
            )

        return self.__class__(
            atom27_xyz=np.concatenate(
                [self.atom27_xyz, np.zeros((n_add, n_atoms, 3))], axis=0
            ),
            atom27_mask=np.concatenate(
                [self.atom27_mask, np.zeros((n_add, n_atoms))], axis=0
            ),
            aatype=np.concatenate(
                [
                    self.aatype,
                    np.full((n_add,), fill_value=len(rc.restype_order), dtype=np.int32),
                ],
                axis=0,
            ),
            b_factors=np.concatenate(
                [self.b_factors, np.zeros((n_add, n_atoms))], axis=0
            ),
            residue_index=np.concatenate(
                [self.residue_index, np.full((n_add,), fill_value=-1, dtype=np.int32)],
                axis=0,
            ),
            chain_index=np.concatenate(
                [self.chain_index, np.zeros((n_add,), dtype=np.int32)], axis=0
            ),
        )