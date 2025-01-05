import os
import matplotlib.pyplot as plt

# Function to plot interaction profiles for base pairs
def plot_interaction_profiles(scores_dict, bins):
    """
    Plot interaction profiles for individual and combined base pair scores.

    Args:
        scores_dict (dict): Dictionary containing the interaction scores for each base pair.
                            Keys are base pair names (e.g., 'AU', 'CG'), and values are lists of scores.
        bins (array-like): Array of bin edges used to categorize distances between base pairs.
                           Length of bins is one more than the length of scores.

    Saves:
        - Individual plots for each base pair in the output directory.
        - A combined plot showing the interaction profiles of all base pairs.
    """
    output_dir = "rna_pdb_files/output_results"
    os.makedirs(output_dir, exist_ok=True)

    # Plot individual interaction profiles for each base pair
    for pair, scores in scores_dict.items():
        plt.figure()
        plt.plot(bins[:-1], scores, marker='o')
        plt.title(f'Interaction Profile for Base Pair {pair}')
        plt.xlabel('Distance (Å)')
        plt.ylabel('Score')
        plt.grid()
        plt.savefig(os.path.join(output_dir, f'{pair}_interaction_profile.png'))
        plt.close()

    # Plot combined interaction profiles for all base pairs
    plt.figure()
    for pair, scores in scores_dict.items():
        plt.plot(bins[:-1], scores, marker='o', label=pair) 

    plt.title('Combined Interaction Profiles')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Score')
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_dir, 'combined_interaction_profile.png'))
    plt.close()
