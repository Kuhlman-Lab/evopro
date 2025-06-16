from Bio.PDB import MMCIFParser, PDBIO
from io import StringIO
import os
import argparse

def cif_to_pdb(mmcif_str):
    parser = MMCIFParser()
    cif_fh = StringIO(mmcif_str) 
    structure = parser.get_structure("structure", cif_fh)
    
    # Truncate residue names to 3 letters
    for model in structure:
        for chain in model:
            for residue in chain:
                if len(residue.resname) > 3:
                    residue.resname = residue.resname[:3]
    
    io = PDBIO()
    io.set_structure(structure)
    output = StringIO()
    io.save(output)
    pdb_string = output.getvalue()
    
    return pdb_string

def convert_cifs_to_pdb(input_dir, output_dir = "pdbs_out"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for filename in os.listdir(input_dir):
        # loop over all .cif files
        if filename.endswith(".cif"):
            input_path = os.path.join(input_dir, filename)
            dest_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.pdb")
            
            with open(input_path, "r") as f:
                content = f.read()
                #print_lines(content)

                # convert to pdb string
            pdb_str = cif_to_pdb(content)
            #print_lines(pdb_str)

            with open(dest_path, "w") as f:
                f.write(pdb_str)
            
            print(f"Saved as: {dest_path}")

    return

def print_lines(str, n=10):
    # prints first n lines of input string
    lines = str.split("\n") # Split the string into a list of lines
    
    for i in range(min(n, len(lines))):
        print(lines[i]) 
    return 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str)
    parser.add_argument("-o", "--output_dir", type=str, default = "pdbs_out")

    args = parser.parse_args()
    
    convert_cifs_to_pdb(**vars(args))