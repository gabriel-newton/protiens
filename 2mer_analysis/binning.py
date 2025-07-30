import pandas as pd
import numpy as np
from tqdm import tqdm

def create_and_save_density_grid(
    input_csv: str, 
    output_npy="angle_density_grid.npy", 
    bins=360):
    """
    Processes a large CSV of angles to create and save a 2D density grid.
    """
    print(f"Creating {bins}x{bins} density grid from '{input_csv}'...")
    
    hist2d = np.zeros((bins, bins), dtype=np.int32)
    
    try:
        # Note: Getting total rows for a remote file can be slow.
        # This is optional and only for the progress bar's total.
        # total_rows = sum(1 for row in open(input_csv)) - 1
        with pd.read_csv(input_csv, chunksize=1_000_000) as reader:
            for chunk in tqdm(reader, desc="Processing Chunks"):
                avg_phi = (chunk['phi_i'] + chunk['phi_i+1']) / 2
                avg_psi = (chunk['psi_i'] + chunk['psi_i+1']) / 2
                
                phi_indices = np.clip(np.round(avg_phi).astype(int) + 180, 0, bins - 1)
                psi_indices = np.clip(np.round(avg_psi).astype(int) + 180, 0, bins - 1)
                
                np.add.at(hist2d, (psi_indices, phi_indices), 1)

    except FileNotFoundError:
        print(f"Error: The file '{input_csv}' was not found.")
        return
        
    np.save(output_npy, hist2d)
    print(f"\nProcessing complete. Density grid saved to '{output_npy}'.")

if __name__ == "__main__":
    # --- UPDATED: Set the full UNC path to your data ---
    source_csv_path = "\\\\rfs01\\rdm01\\GeometryOfProteins\\GabrielNewton\\2mer_angles.csv"
    create_and_save_density_grid(input_csv=source_csv_path)