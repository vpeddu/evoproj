import os
import sys
import csv
import subprocess
import numpy as np
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def convert_bigbed_to_bed(bigbed_dir, output_bed):
    """
    Detects all bigBed files in a directory, converts them to BED format, 
    and merges them into a single BED file.

    Args:
        bigbed_dir (str): Directory containing bigBed files.
        output_bed (str): Path to the output merged BED file.

    Returns:
        bool: True if successful, False otherwise.
    """
    # Find all bigBed files in the directory
    bigbed_files = []
    for file in os.listdir(bigbed_dir):
        if file.endswith('.bb'):
            bigbed_files.append(os.path.join(bigbed_dir, file))
    
    if not bigbed_files:
        print(f"No bigBed files found in {bigbed_dir}")
        return False

    print(f"Found {len(bigbed_files)} bigBed files to process")

    # Process each bigBed file with progress bar
    temp_bed_files = []
    for bigbed in tqdm(bigbed_files, desc="Converting bigBed files"):
        temp_bed = f"{os.path.basename(bigbed)}.bed"
        temp_bed_files.append(temp_bed)
        # Convert bigBed to BED using bigBedToBed
        subprocess.run(["/usr/local/bin/bigBedToBed", bigbed, temp_bed], check=True)

    # Merge all BED files into a single file
    print("Merging BED files...")
    with open(output_bed, 'w') as outfile:
        for temp_bed in temp_bed_files:
            with open(temp_bed, 'r') as infile:
                outfile.write(infile.read())
            os.remove(temp_bed)  # Clean up temporary files
    
    return True

def read_genome_fai(fai_file):
    """
    Reads the genome index file (.fai) and returns a dictionary of chromosome lengths.

    Args:
        fai_file (str): Path to the genome index file.

    Returns:
        dict: A dictionary with chromosome names as keys and their lengths as values.
    """
    chrom_lengths = {}
    with open(fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]               # Chromosome name
            length = int(fields[1])         # Chromosome length
            chrom_lengths[chrom] = length
    return chrom_lengths

def read_bed_file(bed_file, chroms_in_genome):
    """
    Reads a BED file and returns intervals per chromosome.

    Args:
        bed_file (str): Path to the BED file.
        chroms_in_genome (dict): Dictionary of chromosomes present in the genome.

    Returns:
        dict: A dictionary where keys are chromosomes and values are lists of intervals.
    """
    intervals = {}
    print(f"Reading BED file: {bed_file}")
    line_count = 0
    with open(bed_file, 'r') as f:
        for line in f:
            line_count += 1
            if line.startswith(('#', 'track', 'browser')):
                continue  # Skip header lines
            fields = line.strip().split()
            chrom = fields[0]
            if chrom not in chroms_in_genome:
                continue  # Ignore chromosomes not in the genome
            start = int(fields[1])
            end = int(fields[2])
            if chrom not in intervals:
                intervals[chrom] = []
            intervals[chrom].append([start, end])  # Use list for mutability in merging
    print(f"Processed {line_count} lines from BED file")
    return intervals

def sort_and_merge_intervals(intervals):
    """
    Sorts and merges overlapping intervals for each chromosome.

    Args:
        intervals (dict): Dictionary of intervals per chromosome.

    Returns:
        dict: Dictionary with merged intervals per chromosome.
    """
    merged_intervals = {}
    print("Sorting and merging intervals...")
    for chrom in tqdm(intervals.keys(), desc="Processing chromosomes"):
        intervals_chrom = intervals[chrom]
        intervals_chrom.sort(key=lambda x: x[0])  # Sort intervals by start position
        merged = [intervals_chrom[0]]  # Initialize merged list with the first interval
        for start, end in intervals_chrom[1:]:
            last = merged[-1]  # Get the last interval in the merged list
            if start <= last[1]:  # Check if intervals overlap
                last[1] = max(last[1], end)  # Merge intervals
            else:
                merged.append([start, end])  # Add new interval to merged list
        merged_intervals[chrom] = merged
    return merged_intervals

def compute_overlap_bases(intervals_a, intervals_b):
    """
    Computes the total number of overlapping bases between two lists of intervals.

    Args:
        intervals_a (list): List of intervals from Set A for a chromosome.
        intervals_b (list): List of intervals from Set B for the same chromosome.

    Returns:
        int: Total number of overlapping bases.
    """
    total_overlap = 0
    i = j = 0  # Pointers for intervals_a and intervals_b
    len_a = len(intervals_a)
    len_b = len(intervals_b)

    while i < len_a and j < len_b:
        a_start, a_end = intervals_a[i]
        b_start, b_end = intervals_b[j]

        if a_end <= b_start:
            i += 1  # Move to next interval in intervals_a
        elif b_end <= a_start:
            j += 1  # Move to next interval in intervals_b
        else:
            # Calculate overlap
            overlap_start = max(a_start, b_start)
            overlap_end = min(a_end, b_end)
            total_overlap += overlap_end - overlap_start
            # Move pointers
            if a_end <= b_end:
                i += 1
            else:
                j += 1
    return total_overlap

def worker_multi_permutation(args):
    """
    Worker function for multiprocessing that performs multiple permutations.

    Args:
        args (tuple): Contains intervals, chromosome lengths, and the number of permutations.

    Returns:
        list: List of overlap counts for each permutation.
    """
    intervals_a, intervals_b, chrom_lengths, chrom_list, n_perms, worker_id = args
    counts = []  # Store overlap counts
    
    # Create progress bar for this worker
    progress_bar = tqdm(total=n_perms, desc=f"Worker {worker_id}", position=worker_id)
    
    for _ in range(n_perms):
        total_overlap = 0
        for chrom in chrom_list:
            if chrom in intervals_a and chrom in intervals_b:
                chrom_length = chrom_lengths[chrom]
                intervals_a_chrom = intervals_a[chrom]
                # Randomize intervals
                lengths = np.array([end - start for start, end in intervals_a_chrom])
                max_starts = chrom_length - lengths
                random_starts = np.random.randint(0, max_starts + 1, size=len(lengths))
                randomized_intervals_a = np.column_stack((random_starts, random_starts + lengths)).tolist()
                # Merge randomized intervals
                randomized_intervals_a = merge_intervals(randomized_intervals_a)
                # Compute overlap with intervals_b[chrom]
                total_overlap += compute_overlap_bases(randomized_intervals_a, intervals_b[chrom])
        counts.append(total_overlap)  # Append overlap count for this permutation
        progress_bar.update(1)
    
    progress_bar.close()
    return counts

def merge_intervals(intervals):
    """
    Merges overlapping intervals in a list.

    Args:
        intervals (list): List of intervals.

    Returns:
        list: List of merged intervals.
    """
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])  # Ensure intervals are sorted
    merged = [intervals[0]]  # Initialize merged list with the first interval
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1][1] = max(last_end, end)  # Merge intervals
        else:
            merged.append([start, end])  # Add new interval to merged list
    return merged
def write_results_to_csv(output_file, species, repeat, observed_overlap, p_value, n_permutations):
    """
    Writes the permutation test results to a CSV file.
    
    Args:
        output_file (str): Path to the output CSV file.
        species (str): Species name.
        repeat (str): Repeat name.
        observed_overlap (int): Number of overlapping bases observed.
        p_value (str): p-value of the permutation test.
        n_permutations (int): Number of permutations performed.
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Species', 'Repeat', 'Overlap Count', 'p-value', 'Number of Permutations']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerow({
            'Species': species,
            'Repeat': repeat,
            'Overlap Count': observed_overlap,
            'p-value': p_value,
            'Number of Permutations': n_permutations
        })

def main():
    args = sys.argv[1:]
    if len(args) < 5:
        print("Usage: python permutation_test.py path/to/SetA.bed path/to/SetB_bigbed_directory path/to/genome.fa.fai species repeat [n_permutations]")
        sys.exit(1)

    # Parse command-line arguments
    bed_a_file = args[0]
    bigbed_dir = args[1]
    fai_file = args[2]
    species = args[3]
    repeat = args[4]
    tf = args[5]
    n_permutations = int(args[6]) if len(args) > 5 else 10000

    print(f"Starting permutation test with {n_permutations} permutations")
    print(f"Species: {species}, Repeat: {repeat}")
    
    # Convert bigBed files to a single merged BED file
    merged_bed_b_file = "merged_SetB.bed"
    if not convert_bigbed_to_bed(bigbed_dir, merged_bed_b_file):
        print("Error processing bigBed files. Exiting.")
        sys.exit(1)

    # Read chromosome lengths from genome index file
    print("Reading genome index file...")
    chrom_lengths = read_genome_fai(fai_file)
    
    # Read intervals from BED files
    intervals_a = read_bed_file(bed_a_file, chrom_lengths)
    intervals_b = read_bed_file(merged_bed_b_file, chrom_lengths)

    # Sort and merge intervals
    intervals_a = sort_and_merge_intervals(intervals_a)
    intervals_b = sort_and_merge_intervals(intervals_b)

    # Compute observed overlap
    print("Computing observed overlap...")
    observed_overlap = 0
    for chrom in tqdm(chrom_lengths.keys(), desc="Processing chromosomes"):
        if chrom in intervals_a and chrom in intervals_b:
            observed_overlap += compute_overlap_bases(intervals_a[chrom], intervals_b[chrom])
    
    print(f"Observed overlap: {observed_overlap} bases")
    print(f"Starting permutation testing with {n_permutations} permutations...")

    # Set up multiprocessing
    n_processes = min(cpu_count(), n_permutations)  # Limit to 8 processes max for display purposes
    print(f"Using {n_processes} parallel processes")
    pool = Pool(processes=n_processes)

    # Distribute permutations evenly across processes
    perms_per_process = n_permutations // n_processes
    extra = n_permutations % n_processes
    # Create a list with the number of permutations each process should perform
    perms_distribution = [perms_per_process + (1 if i < extra else 0) for i in range(n_processes)]

    chrom_list = list(chrom_lengths.keys())  # List of chromosomes
    # Prepare arguments for each process
    args_list = [(intervals_a, intervals_b, chrom_lengths, chrom_list, n_perms, i) 
                 for i, n_perms in enumerate(perms_distribution)]

    # Run permutations in parallel
    results = pool.map(worker_multi_permutation, args_list)
    pool.close()
    pool.join()

    # Flatten results from all processes into a single list
    print("Aggregating results...")
    overlaps = [overlap for sublist in results for overlap in sublist]

    # Calculate p-value
    greater_or_equal = sum(1 for overlap in overlaps if overlap >= observed_overlap)
    p_value = greater_or_equal / n_permutations

    # Adjust p-value formatting
    if p_value == 0.0 and greater_or_equal == 0:
        p_value_str = f"< {1 / n_permutations:.16f}"
    else:
        p_value_str = f"{p_value:.16f}"

    # Print the results in the required format
    print("\n" + "=" * 50)
    print(f"Number of overlapping bases observed: {observed_overlap}")
    print(f"p-value: {p_value_str}")
    print("=" * 50)

    # Write results to CSV
    output_csv = f"{species}_{repeat}_{tf}_permutation_results.csv"
    write_results_to_csv(output_csv, species, repeat, observed_overlap, p_value_str, n_permutations)
    print(f"Results written to {output_csv}")

if __name__ == "__main__":
    main()