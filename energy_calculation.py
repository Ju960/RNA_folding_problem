from Training_script import extract_c3_prime_coordinates

# Function to interpolate the score for a given distance and base pair
def interpolate_score(dist, base_pair):
    """
    Interpolate the interaction score for a given distance and base pair.

    Args:
        dist (float): The distance between two C3' atoms in Ångströms.
        base_pair (str): The base pair (e.g., 'AU', 'CG') for which the score is required.

    Returns:
        float: Interpolated score for the given distance and base pair.
               Returns None if no score file is found or the distance is outside the range.
    """
    score_file = f'rna_pdb_files/output_results/{base_pair}_log_ratio.txt'
    intervals = []
    scores = []   

    try:
        # Read the score file and extract intervals and scores
        with open(score_file, 'r') as f:
            for line in f:
                distance_range, score = line.strip().split(':')
                distance_range = distance_range.replace('Å', '').strip()
                min_dist, max_dist = map(float, distance_range.split('-'))
                score = float(score.strip())
                intervals.append((min_dist, max_dist))
                scores.append(score)
    except FileNotFoundError:
        print(f"Error: Score file not found for {base_pair}")
        return None

    # Find the interval to which the distance belongs
    for i in range(len(intervals)):
        min_dist, max_dist = intervals[i]
        if min_dist <= dist <= max_dist:  # Check if the distance is within the current interval
            score1 = scores[i]
            # Check if interpolation with the next interval is possible
            if i + 1 < len(intervals):
                min_dist_next, max_dist_next = intervals[i + 1]
                score2 = scores[i + 1]
                
                # Perform linear interpolation between the current and next intervals
                interpolated_score = score1 + (dist - min_dist) / (max_dist - min_dist) * (score2 - score1)
                return interpolated_score
            else:
                return score1

    return None


# Function to extract distances between C3' atoms and associate them with base pairs
def extract_c3_prime_distances_with_pairs(pdb_file, output_file):
    """
    Extract distances between C3' atoms and associate them with their corresponding base pairs.

    Args:
        pdb_file (str): Path to the input PDB file containing RNA structure.
        output_file (str): Path to the output file to save the distances, base pairs, and scores.

    Returns:
        float: Total Gibbs free energy calculated as the sum of scores for all base pairs.
    """
    # Call the function to extract C3' coordinates from the PDB file
    base_pairs_coords, all_pairs_coords = extract_c3_prime_coordinates(pdb_file)
    
    distances_with_pairs = []
    total_energy = 0.0

    # Calculate distances and associate them with base pairs
    for dist, pair in all_pairs_coords:
        # Interpolate the score based on the distance and base pair
        score = interpolate_score(dist, pair)
            
        if score is not None:
            total_energy += score
            distances_with_pairs.append((dist, pair, score))
        else:
            distances_with_pairs.append((dist, pair, "N/A"))

    # Write the results to the output file
    with open(output_file, 'w') as out_f:
        out_f.write("Distance (Å)\tBase Pair\tScore\n")
        for dist, pair, score in distances_with_pairs:
            if score != "N/A":
                out_f.write(f"{dist:.3f}\t{pair}\t{score:.3f}\n")
            else:
                out_f.write(f"{dist:.3f}\t{pair}\t{score}\n")

    return total_energy
