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

PLEASE MAKE SURE TO USE "STABLE" BRANCH for code used in preprint. Current working version is on "dev branch".

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

## HOW TO WRITE AN EVOPRO SCORING FUNCTION
The score functions are a way to generate a “score” from the AlphaFold2 result for each sequence of the optimization pool – and this score is used for ranking and filtering the sequences (so it is only the relative value of the score that matters). Each score function therefore takes in the “results” dictionary and parses the different components such as the predicted structure in PDB form, pLDDT, PAE, etc. The score function should then return a final overall score for the sequence that can contain different weighted score components.
 
EvoPro includes a few functions that are prewritten to be called within your score function, that calculate different score metrics which may then be used as weighted components in the final overall score that is returned from your score function. For example, a function to calculate the average pLDDT over a predicted PDB (score_plddt_confidence) or to calculate the contacts at an interface, where each contact is weighted by the PAE of the interaction (score_contacts_pae_weighted). You can import these from evopro.score_funcs.score_funcs, as is done at the top of score_binder.py, or use your own functions instead.
Important to know: the more negative scores are better, so we negate score terms that we want to maximize and keep positive the score terms that we want to minimize.
 
See below for an example score function with annotated code: 
 
def score_binder(results, dsobj, contacts=None, distance_cutoffs=None):
#This is the primary scoring function that is called in the main code. The purpose of this
#function is to determine whether it is the complex prediction (with both binder and
#target) or the monomeric (just binder) prediction that is being evaluated and return the
#corresponding score. The main code will then pair the complex and monomer predictions and
#sum the scores from each.
 
    from alphafold.common import protein
    #Generate the pdb from the alphafold output
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    #Parse chains and residues from the pdb
    #chains will look like [“A”, “B”] and residues like [“A1”, “A2”, … , “B1”, “B2”, …]
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        #If the prediction has more than one chain call the complex scoring function.
        #Here, contacts and distance_cutoffs come as arguments from the main code
        #(see more descriptions of both inside this function below)
        return score_binder_complex(results, dsobj, contacts, distance_cutoffs)
    else:
        #If the prediction has only one chain call the monomer scoring function.
        return score_binder_monomer(results, dsobj)
 
def score_binder_complex(results, dsobj, contacts, distance_cutoffs):
#This is the complex scoring function that is called from the primary score function.
 
#Here, contacts is a tuple containing three items which all come from user inputs:
#(interface_contacts, bonus_contacts, penalty_contacts)
#Interface contacts come from the flag --define_contact_area where the user defines
#which residues on the target protein are being targeted for binding.
#Bonus contacts come from the flag -- bonus_contacts (optional) where the user defines
#residues that are given a “bonus” if contacted by the binder.
#For this, we usually use a single residue in the middle of the target interface to help
#select for binders that are in the right place.
#Penalty contacts come from the flag -- penalize_contacts (optional) where the user defines
#residues that are given a “penalty” if contacted by the binder.
 
#Distance cutoffs also contains 3 elements representing the Angstrom cutoff for each
#of these types of contacts. The default is (4,4,8) and can be changed through the flags.
 
 
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
 
    #setting defaults if no user input
    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
 
    #Get interface contacts on target protein
    reslist1 = contacts[0]
    #Get residues on binder (here, all of chain B) that can bind to the interface contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
   
    #Get list of interacting residues and PAE-weighted contact score
    contact_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0])
   
    bonuses = 0
    #Get bonus contacts on target protein
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, bonus_contactscore = score_contacts_pae_weighted(results, pdb, bonus_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1])
        for contact in bonus_contacts:
            if contact[0][0:1] == 'A' and int(contact[0][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[0]))
            if contact[1][0:1] == 'A' and int(contact[1][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[1]))
       
    #Give each bonus contact a score value of -3.
    bonus = -bonuses * 3
   
    penalties = 0
    penalty_resids = contacts[2]
    if penalty_resids:
        penalty_contacts, penalty_contactscore = score_contacts_pae_weighted(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2])
        for contact in penalty_contacts:
            if contact[0][0:1] == 'A' and int(contact[0][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[0]))
            if contact[1][0:1] == 'A' and int(contact[1][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[1]))
     
    #Give each penalty contact a score value of +3
    penalty = penalties * 3
   
    #Get the average PAE of interface contacts (pae_interaction)
    num_contacts = len(contacts)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
   
    #Calculate the weighted score to return
    score = -contactscore + penalty + bonus
 
    #Return values. Here, first term is the overall score to use for ranking and second
    #term has the individual score components that will be printed out in the log file.
    #NOTE: pdb, results MUST be returned and should always be the last 2 elements.
 
    return score, (score, len(contacts), contactscore, pae_per_contact, bonus, penalty), contacts, pdb, results
 
def score_binder_monomer(results, dsobj):
#This is the monomer scoring function that is called from the primary score function.
 
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
 
    #Parse chains and residues from the pdb
    chains, residues, resindices = get_coordinates_pdb(pdb)
 
    #Select residues to be included in average pLDDT score
    #Here, we select all residues in the binder-only AF2 prediction
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
 
    #calculate overall score to be returned from this function
    #We divide by 10 since the values are in generally the range of 80-100
    #so this term does not dominate over the complex scores (in the 10s range)
    score = -confscore2/10
 
    return score, (score, confscore2), pdb, results
 

<!-- CONTACT -->
## Contact

Amrita Nallathambi - [@amritanalla](https://twitter.com/twitter_handle) - anallathambi@unc.edu

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[conda-install-link]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
