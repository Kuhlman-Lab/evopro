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

EvoPro is a genetic algorithm-based protein binder optimization pipeline, used for "in silico evolution" of highly accurate, tight protein binders.

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

4. Set paths at the top of the script evopro/run/run_geneticalg_gpus.py.
   

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Usage of EvoPro requires GPUs.
The running directory should contain input files for the target protein, starting scaffolds for binders, and flags files for EvoPro specifications and AF2 specifications. See evopro/examples/ for examples of directory setups.

The residue_specs.json file should be generated from pdb or sequence file.
Specify pdb/sequence file and residues to mutate in json.flags.
Then, in the running directory:
```sh
python /path/to/evopro/run/generate_json.py @json.flags
 ```


You can set options for EvoPro in the evopro.flags file.
 To run EvoPro, confirm availability of GPUs. Then, 
 ```sh
conda activate evopro
python /path/to/evopro/run/run_geneticalg_gpus.py @evopro.flags
 ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Amrita Nallathambi - [@amritanalla](https://twitter.com/twitter_handle) - anallathambi@unc.edu

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[conda-install-link]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
