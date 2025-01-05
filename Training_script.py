import os
import numpy as np
import math

def extract_c3_prime_coordinates(pdb_file):
    """
    Extracts the coordinates of C3' atoms from a PDB file and calculates distances between them
    to form base pairs that respect the i and i+4 separation constraint.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        tuple: 
            - base_pairs (dict): Dictionary containing lists of distances for each base pair (e.g., 'AU', 'GC').
            - all_pairs (list): List of tuples containing distances and corresponding base pairs.
    """
    base_pairs = {
        'AA': [], 'AU': [], 'AC': [], 'AG': [], 'UU': [], 'CU': [], 'GU': [], 'CC': [], 'CG': [], 'GG': []
    }
    all_pairs = []  # List of tuples (distance, base pair)
    c3_prime_coords = []  # List of tuples (residue name, coordinates)

    # Read the PDB file line by line
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            # Check if the line describes an atom or heteroatom
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                if atom_name == "C3'":  # Only consider C3' atoms
                    chain = line[21]
                    residue_id = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    c3_prime_coords.append((residue_name, (x, y, z)))

    # Compute distances between C3' atoms that are separated by at least 4 residues
    for i in range(len(c3_prime_coords)):
        base1, coord1 = c3_prime_coords[i]
        for j in range(i + 4, len(c3_prime_coords)):
            base2, coord2 = c3_prime_coords[j]
            distance = np.linalg.norm(np.array(coord1) - np.array(coord2))
            base_pair = ''.join(sorted([base1, base2]))  # Create a valid base pair name

            if base_pair in base_pairs:
                base_pairs[base_pair].append(distance)
            all_pairs.append((distance, base_pair))

    return base_pairs, all_pairs


def calculate_observed_frequencies(all_base_pairs):
    """
    Calculate the observed frequencies of base pair distances and the reference counts
    across specified distance bins.

    Args:
        all_base_pairs (dict): Dictionary containing lists of distances for each base pair.

    Returns:
        tuple:
            - observed_frequencies (dict): Observed frequency counts for each base pair.
            - reference_counts (np.ndarray): Total counts of all distances across bins.
            - bins (np.ndarray): Bin edges used for the histogram.
    """
    bins = np.linspace(0, 20, 21)  # Define distance bins from 0 to 20 Å with a step of 1 Å
    observed_frequencies = {pair: np.zeros(len(bins) - 1) for pair in all_base_pairs.keys()}
    reference_counts = np.zeros(len(bins) - 1)

    # Calculate frequency counts for each base pair and all distances combined
    for pair, distances in all_base_pairs.items():
        counts, _ = np.histogram(distances, bins=bins)
        observed_frequencies[pair] += counts

    all_distances = [d for distances in all_base_pairs.values() for d in distances]
    reference_counts, _ = np.histogram(all_distances, bins=bins)

    return observed_frequencies, reference_counts, bins


def compute_log_ratio(observed_frequencies, reference_counts, bins, max_score=10):
    """
    Compute the log-ratio between observed and reference frequencies for each base pair
    and save the results to output files.

    Args:
        observed_frequencies (dict): Observed frequency counts for each base pair.
        reference_counts (np.ndarray): Total counts of all distances across bins.
        bins (np.ndarray): Bin edges used for the histogram.
        max_score (float): Maximum score value to cap the log-ratio.

    Returns:
        dict: Log-ratio scores for each base pair.
    """
    scores_dict = {}
    output_dir = "rna_pdb_files/output_results"
    os.makedirs(output_dir, exist_ok=True)

    # Calculate reference frequency distribution
    total_reference = np.sum(reference_counts)
    reference_frequency = reference_counts / total_reference if total_reference > 0 else np.zeros_like(reference_counts)

    # Compute log-ratio for each base pair
    for pair, observed_counts in observed_frequencies.items():
        observed_total = np.sum(observed_counts)
        observed_frequency = observed_counts / observed_total if observed_total > 0 else np.zeros_like(observed_counts)

        scores = []
        for obs, ref in zip(observed_frequency, reference_frequency):
            if ref > 0 and obs > 0:
                score = -math.log(obs / ref)
                score = min(score, max_score)
            else:
                score = max_score
            scores.append(score)

        scores_dict[pair] = scores

        # Save the scores to a file
        filename = os.path.join(output_dir, f"{pair}_log_ratio.txt")
        with open(filename, 'w') as f:
            for i in range(len(scores)):
                f.write(f"{bins[i]:.1f}-{bins[i+1]:.1f} Å : {scores[i]:.3f}\n")
        print(f"File generated: {os.path.abspath(filename)}")

    return scores_dict


def process_pdb_files(input_dir):
    """
    Process all PDB files in the input directory, extract distances between base pairs,
    and compute the log-ratio for each pair.

    Args:
        input_dir (str): Path to the directory containing PDB files.

    Returns:
        tuple:
            - scores_dict (dict): Log-ratio scores for each base pair.
            - bins (np.ndarray): Bin edges used for the histogram.
    """
    all_base_pairs = {pair: [] for pair in ['AA', 'AU', 'AC', 'AG', 'UU', 'CU', 'GU', 'CC', 'CG', 'GG']}
    
    # Process each PDB file in the directory
    for pdb_filename in os.listdir(input_dir):
        if pdb_filename.endswith(".pdb"):
            pdb_file = os.path.join(input_dir, pdb_filename)
            print(f"Processing file {pdb_file}...")

            # Extract base pair distances
            base_pairs_coords, all_pairs_coords = extract_c3_prime_coordinates(pdb_file)

            # Append distances to the global list for each base pair
            for pair in all_base_pairs.keys():
                all_base_pairs[pair].extend(base_pairs_coords[pair])

    # Calculate observed and reference frequencies
    observed_frequencies, reference_counts, bins = calculate_observed_frequencies(all_base_pairs)

    # Compute log-ratio and save results
    scores_dict = compute_log_ratio(observed_frequencies, reference_counts, bins)
    return scores_dict, bins
