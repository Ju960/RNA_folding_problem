# RNA_FOLDING_PROBLEM

`RNA_folding_problem` is a Python package designed for the RNA structures analyses, focusing on calculating the interaction profiles between base pairs in RNA. The analysis includes calculating the distances between C3' atoms of ribonucleotides and computing log-ratio scores for each base pair using PDB files. The results can be visualized and saved for further analysis. 


## Package structure

RNA_INTERACTION_ANALYSIS/  
├── [__init__.py]                    # Initialization script for the package  
├── [rna_pdb_files]                  # Folder containing PDB files and output results  
    ├── [output_results]             # Folder for saving computed results (e.g., interaction profiles, scores)  
├── [Native_Pred_structure]          # Folder containing native and predictive PDB files with their scores  
├── [Training_script.py]             # Functions to extract C3' coordinates from PDB files, compute distances, frequencies and log-ratio scores for each base pair  
├── [Plot_interaction_score.py]      # Function to plot interaction profiles and combined results  
├── [energy_calculation.py]          # Functions to calculate the Gibbs Energy for 1 native and 1 predictive  
├── [README.md]                      # Documentation for the package  
├── [requirements.txt]               # Required Python packages for the project  
└──[example.ipynb]                   # Jupyter notebook example for using the package  


## Installation

Clone this repository and install the required dependencies:

```bash
git clone <https://github.com/Ju960/RNA_folding_problem.git>
cd RNA_folding_problem
pip install -r requirements.txt
```

## Usage

### Extracting C3' Coordinates and Calculating Distances

To extract the C3' atom coordinates from a folder containing PDB files, calculate distances, frequencies and scores for each base pairs:

```python
from RNA_folding_problem.Training_script import process_pdb_files

pdb_file = "path/to/your/pdb_folder"
scores_dict, bins = process_pdb_files(input_dir)

# scores_dict is Log-ratio scores for each base pair,
# bins contains the bins from 0 to 20 Å in step 1 Å.

```

### Plotting Interaction Profiles

To plot and save the interaction profiles for each base pair:

```python
from RNA_folding_problem.Plot_interaction_score import plot_interaction_profiles
plot_interaction_profiles(scores_dict, bins)

```

### Calculating Gibbs Energy

To calculate the Gibbs Energy for a given PDB file:

```python
from RNA_folding_problem.energy_calculation import extract_c3_prime_distances_with_pairs

pdb_file = "path/to/your/file.pdb"
output_file = "path/to/output_file.txt"
total_energy = extract_c3_prime_distances_with_pairs(pdb_file, output_file)

print(f"The total Gibbs free energy is: {total_energy}")

```

## Examples

Detailed usage examples can be found in the [example.ipynb](../example.ipynb) file

## Authors

Julia GOUNIN - M2 GENIOMHE University of Evry

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
