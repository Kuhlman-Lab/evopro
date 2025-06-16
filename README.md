<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

<!-- ABOUT THE PROJECT -->
## About EvoPro

EvoPro is a genetic algorithm-based protein binder optimization pipeline, used in published work for in silico evolution of highly accurate, tight protein binders.  
Now with updated code including multistate design - PREPRINT COMING SOON!

PLEASE MAKE SURE TO USE "MAIN" BRANCH for stable version. Branch "pd1_paper" contains old code used in paper. Current working version is on "dev branch". 

[paper]  


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

Steps to set up a copy of EvoPro locally are detailed below.

### Prerequisites

Installation of Anaconda is required to load dependencies.
* Installing conda: [conda-install-link]

### EvoPro Installation

1. Clone the repo:
   ```sh
   git clone https://github.com/Kuhlman-Lab/evopro.git 
   ```
2. Clone our AF2 and ProteinMPNN repos:
   ```sh
   git clone https://github.com/Kuhlman-Lab/alphafold.git
   git clone https://github.com/Kuhlman-Lab/proteinmpnn.git
   ```
3. Load AlphaFold2 model weights from source using script: https://github.com/Kuhlman-Lab/alphafold/blob/main/setup/download_alphafold_params.sh 

4. Set up conda environment:
   ```sh
   conda env create -n evopro -f setup_conda.yaml
   pip3 install --upgrade jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
   python3 -m pip install /path/to/alphafold/alphafold/
   ```

5. Set the path to your EvoPro installation at the top of run/generate_json.py and run/run_evopro.py.


<!-- USAGE EXAMPLES -->
## Usage

Usage of EvoPro requires GPUs.
The running directory should contain a sequence file with input sequences, yaml file for EvoPro specifications, and flags files for AF2 and ProteinMPNN specifications. See below for flag options and evopro/examples/ for examples of directory setups. 

The residue_specs.json file should be generated from the sequence file (can also be generated from a PDB).
Specify sequence file and which residues to mutate in json.flags. (Include symmetry here if needed)
Then, in the running directory:
```sh
python /path/to/evopro/run/generate_json.py @json.flags
 ```

You can set options for EvoPro in the evopro.yaml file.
 To run EvoPro, confirm availability of GPUs. Then, 
 ```sh
conda activate evopro
python /path/to/evopro/run/run_evopro.py --config_file evopro.yaml
 ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Generating the JSON file

Specify your flag options in json.flags, and provide either a PDB or sequence file to extract starting sequences from. All flag options are listed below:  
--pdb  
default=None, type=str,  
Path to and name of PDB file to extract chains and sequences. Alternatively use sequence file (below).

--sequence_file  
default=None, type=str,  
Path to and name of text file to extract chains and sequences. Only provide if there is no PDB file.

--mut_res  
default=''”, type=str,  
PDB chain and residue numbers to mutate, separated by commas. (Make sure there are no extraneous spaces between commas)  
Here, you can use "\*" to specify the whole chain (eg. A\*), or "<" to specify all residues in a chain after a specific residue (eg B17<).

--default_mutres_setting  
default='all', type=str,  
Default setting for residues to mutate. Individual residues can be changed manually after generation of the file. Default is all. Other examples are “all-PC” (allow any mutation except proline and cysteine), “AVILMFYW” (specify which amino acids are allowed).  

--output  
default='', type=str,  
path and name of output json file

--symmetric_res  
default='', type=str,  
PDB chain and residue numbers to force symmetry separated by a colon

## Specifying the AlphaFold2 options

Specify your AF2 options in af2.flags.  
Notably, we suggest turning off MSA generation or providing pre-computed MSA and templates so that the AF2 runner will not have to query MMseqs server for every prediction (this will make runtime much longer and may cause errors).  

To turn off MSA generation:  
--msa_mode single_sequence  

To use precomputed MSA:  
--msa_mode single_sequence  
--custom_msa_path /path/to/directory/with/a3m/files/

To use custom templates:  
--use_templates  
--custom_template_path templates  

Since we have modified the original AF2 code to allow for custom template databases, you will have to make sure each .pdb file within the templates folder has a name consisting of 4 letters and numbers, essentially mimicking a file from the PDB database with a PDB code (although the file name does not actually have to be a real PDB code). Some examples could be “temp.pdb”, “1uzg.pdb”, “PPPP.pdb”, etc.


<!-- CONTACT -->
## Contact

Amrita Nallathambi - anallathambi@unc.edu  
Since I am a grad student, this code may not be perfect. Please email me with any questions or concerns!

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[conda-install-link]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html  
[preprint]: https://www.biorxiv.org/content/10.1101/2023.05.03.539278v1
[paper]: https://www.pnas.org/doi/10.1073/pnas.2307371120
