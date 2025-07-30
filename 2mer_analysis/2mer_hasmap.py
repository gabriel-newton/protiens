import pandas as pd
from collections import Counter
from tqdm import tqdm

filepath = "\\\\rfs01\\rdm01\\GeometryOfProteins\\GabrielNewton\\2mer_angles.csv"

def find_angle_peaks(filepath=filepath, chunksize=1_000_000, top_n=50):
    """
    Reads a large CSV of 2-mer angles, finds the most frequent
    discretized angle combinations, and prints the top N results.

    Args:
        filepath (str): Path to the input CSV file.
        chunksize (int): Number of rows to read into memory at a time.
        top_n (int): The number of top peaks to display.
    """
    print(f"Starting analysis of '{filepath}'...")
    
    # Initialize a Counter object to store the frequency of each 4-angle tuple.
    angle_counts = Counter()
    
    # Use pandas chunking to iterate through the large file without
    # loading it all into memory. A progress bar is added via tqdm.
    try:
        # Get total rows for tqdm progress bar
        total_rows = sum(1 for row in open(filepath)) - 1 # -1 for header
        
        with pd.read_csv(filepath, chunksize=chunksize) as reader:
            for chunk in tqdm(reader, total=total_rows//chunksize, desc="Processing Chunks"):
                
                # 1. Discretize: Round angles to the nearest integer.
                # .values converts the DataFrame chunk to a NumPy array for speed.
                discretized_chunk = chunk.round().astype(int).values
                
                # 2. Count: Update the master counter with counts from this chunk.
                # We convert each row (a numpy array) into a tuple so it can be a dictionary key.
                angle_counts.update(map(tuple, discretized_chunk))

    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return

    print("\nAnalysis complete. Identifying the highest peaks...")

    # 3. Find Peaks: Get the most frequent angle combinations.
    most_common_angles = angle_counts.most_common(top_n)

    print(f"\n--- Top {top_n} Angle Combination Peaks ---")
    print("      (phi_i, psi_i, phi_i+1, psi_i+1)   |   Count")
    print("---------------------------------------------|-----------")
    for (angles, count) in most_common_angles:
        print(f" {str(angles):<36} | {count}")

# --- Main execution ---
if __name__ == "__main__":
    find_angle_peaks()