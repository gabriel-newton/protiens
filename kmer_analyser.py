import pandas as pd
import os
import glob
import ast
from collections import Counter
from typing import List, Dict, Tuple
from tqdm import tqdm

# This initializes tqdm for use with pandas .progress_apply()
tqdm.pandas(desc="Processing Rows")

# === Global Parameters ===
# Note: These Windows-style UNC paths will fail on Linux unless the remote share
# is mounted. For example, if \\rfs01\rdm01 is mounted at /mnt/rfs01, the path
# should be changed to /mnt/rfs01/GeometryOfProteins/...
GEO_DATA_PATH = "\\\\rfs01\\rdm01\\GeometryOfProteins\\Ziqiu\\invariant_result\\geo_info_18Mar2024.csv"
INVARIANT_DATA_PATH = "\\\\rfs01\\rdm01\\GeometryOfProteins\\Ziqiu\\invariant_result\\invariant_with_weak\\"

class KmerAnalyzer:
    """
    A class to handle the analysis of k-mers from protein geometric sequences.
    """
    def __init__(self, k: int, depth: int, geo_data_path: str = GEO_DATA_PATH, invariant_data_path: str = INVARIANT_DATA_PATH):
        """
        Initializes the analyzer, using global paths as defaults.
        """
        self.k = k
        self.depth = depth
        self.geo_data_path = geo_data_path
        self.invariant_data_path = invariant_data_path
        
        try:
            self.geo_df = pd.read_csv(self.geo_data_path)
            # Magic number -3 removes the trailing 'END' from sequences.
            self.sequences = self.geo_df['geo_seq'].dropna().astype(str).str[:-3].tolist()
        except FileNotFoundError:
            print(f"Error: The file '{self.geo_data_path}' was not found.")
            self.geo_df = pd.DataFrame()
            self.sequences = []
            
        self.df_with_locations = None
        self._file_cache = {}
        
        # --- OPTIMIZATION: Build the file map once during initialization ---
        self.invariant_file_map = {}
        self._create_invariant_file_map()

    def _create_invariant_file_map(self):
        """
        OPTIMIZATION: Scans the invariant data directory once and creates a mapping 
        from a unique identifier tuple to the full file path. This avoids repeated
        file system searches in the main processing loop.
        """
        print("Building invariant file map... this may take a moment.")
        search_path = os.path.join(self.invariant_data_path, "*.csv")
        all_files = glob.glob(search_path)

        for f_path in tqdm(all_files, desc="Mapping invariant files"):
            filename = os.path.basename(f_path)
            parts = filename.split('-')
            if len(parts) >= 4:
                try:
                    pdb_id = parts[0]
                    residue_id = parts[1]
                    model_id = parts[2]
                    chain_id = parts[3]
                    key = (pdb_id, residue_id, model_id, chain_id)
                    
                    if key not in self.invariant_file_map:
                        self.invariant_file_map[key] = f_path
                except IndexError:
                    pass
        print(f"File map built. Found {len(self.invariant_file_map)} unique invariant files.")

    def find_top_kmers(self) -> List[Tuple[str, int]]:
        """
        Finds the most frequent k-mers and their counts across all sequences.
        Returns a list of (k-mer, count) tuples.
        """
        all_kmers = []
        for s in tqdm(self.sequences, desc=f"Finding all {self.k}-mers"):
            if len(s) >= self.k:
                for i in range(len(s) - self.k + 1):
                    all_kmers.append(s[i:i+self.k])
        
        if not all_kmers:
            return []

        frequencies = Counter(all_kmers)
        top_kmers_with_counts = frequencies.most_common(self.depth)
        return top_kmers_with_counts

    def _find_single_kmer_locations(self, sequence: str, kmer_to_find: str) -> List[int]:
        """Helper function to find all locations of a specific k-mer in a sequence."""
        locations = []
        k = len(kmer_to_find)
        if not isinstance(sequence, str):
            return locations
            
        for i in range(len(sequence) - k + 1):
            if sequence[i:i+k] == kmer_to_find:
                locations.append(i)
        return locations

    def create_location_data(self, force_rerun: bool = False):
        """
        Creates and saves a DataFrame with location information for the top k-mers.
        """
        output_dir = f"k{self.k}"
        output_filename = os.path.join(output_dir, f"geo_info_k{self.k}_locations.csv")
        os.makedirs(output_dir, exist_ok=True)

        run_analysis = True
        if os.path.exists(output_filename) and not force_rerun:
            header_df = pd.read_csv(output_filename, nrows=0)
            kmer_columns_in_file = [col for col in header_df.columns if len(col) == self.k]
            
            if len(kmer_columns_in_file) >= self.depth:
                print(f"Loading existing data from '{output_filename}'...")
                self.df_with_locations = pd.read_csv(output_filename)
                run_analysis = False
            else:
                print(f"File only contains data for depth {len(kmer_columns_in_file)}. Re-running for depth {self.depth}...")

        if run_analysis:
            print("Running full analysis to generate location data...")
            top_kmers_with_counts = self.find_top_kmers()
            if not top_kmers_with_counts:
                self.df_with_locations = pd.DataFrame()
                return self.df_with_locations

            top_kmers = [kmer for kmer, count in top_kmers_with_counts]
            self.df_with_locations = self.geo_df.copy()
            
            for kmer in tqdm(top_kmers, desc="Finding k-mer locations"):
                self.df_with_locations[kmer] = self.df_with_locations['geo_seq'].dropna().astype(str).apply(
                    lambda seq: self._find_single_kmer_locations(sequence=seq, kmer_to_find=kmer)
                )
            
            mask = self.df_with_locations[top_kmers].apply(
                lambda row: any(isinstance(loc_list, list) and len(loc_list) > 0 for loc_list in row),
                axis=1
            )
            self.df_with_locations = self.df_with_locations[mask]

            self.df_with_locations.to_csv(output_filename, index=False)
            print(f"Saved new location data to '{output_filename}'.")
        
        return self.df_with_locations

    def _load_single_residue_file(self, row: pd.Series) -> Tuple[pd.DataFrame | None, str | None]:
        """
        MODIFIED: Finds and loads the invariant data file for a row using the 
        pre-built file map, with caching. 
        Returns a tuple of (DataFrame, filepath).
        """
        try:
            key = (
                row['pdb_id'],
                str(row['residue_id']),
                str(row['model_id']),
                row['chain_id']
            )
        except KeyError:
            return None, None

        filepath = self.invariant_file_map.get(key)

        if not filepath:
            return None, None

        if filepath in self._file_cache:
            # Return cached dataframe and the filepath
            return self._file_cache[filepath], filepath

        try:
            data_df = pd.read_csv(filepath)
            self._file_cache[filepath] = data_df 
            return data_df, filepath
        except Exception: 
            self._file_cache[filepath] = None 
            return None, None

    def _save_kmer_data(self, kmer_data_collection: Dict[str, List], output_dir: str):
        """Helper method to save the collected k-mer data to CSV files."""
        for kmer, data_list in kmer_data_collection.items():
            if data_list:
                filepath = os.path.join(output_dir, f"{kmer}.csv")
                write_header = not os.path.exists(filepath)
                df_to_save = pd.DataFrame(data_list)
                df_to_save.to_csv(filepath, mode='a', header=write_header, index=False)

    def extract_invariant_data(self, force_rerun: bool = False, checkpoint_interval: int = 100) -> Dict[str, pd.DataFrame]:
        """
        MODIFIED: Extracts the 9-part invariant data for each top k-mer.
        - occurrence_id is now local to each source file.
        - source_file column now contains the full path to the data file.
        """
        if self.df_with_locations is None:
            print("Location data not found. Running .create_location_data() first...")
            self.create_location_data()
            if self.df_with_locations is None or self.df_with_locations.empty:
                print("Failed to create location data. Aborting.")
                return {}

        output_dir = f"k{self.k}"
        os.makedirs(output_dir, exist_ok=True)
        progress_file = os.path.join(output_dir, '.progress.csv')
        
        df_to_process = self.df_with_locations.copy()
        processed_indices = set()

        if os.path.exists(progress_file) and not force_rerun:
            print(f"Found progress file: '{progress_file}'. Resuming analysis.")
            try:
                progress_df = pd.read_csv(progress_file, header=None, names=['processed_index']).dropna()
                processed_indices = set(progress_df['processed_index'].astype(int))
                df_to_process = df_to_process[~df_to_process.index.isin(processed_indices)]
            except (pd.errors.EmptyDataError, KeyError):
                print("Progress file is empty or malformed. Starting from scratch.")
                if os.path.exists(progress_file):
                    os.remove(progress_file)

        if df_to_process.empty:
            print("All rows have already been processed. Nothing to do.")
            # Still load and return existing data
            return self._load_final_data(output_dir)

        all_top_kmers = [col for col in self.df_with_locations.columns if len(col) == self.k and isinstance(col, str)]
        
        temp_kmer_data_buffer = {kmer: [] for kmer in all_top_kmers}
        
        invariant_cols = [
            'length(N)', 'length(A)', 'length(C)', 
            'angle(N)', 'angle(A)', 'angle(C)',
            'tau(NA)', 'tau(AC)', 'tau(CN)'
        ]

        processed_in_session = []
        
        try:
            for index, row in tqdm(df_to_process.iterrows(), total=len(df_to_process), desc="Processing files"):
                # MODIFICATION: Unpack tuple of (DataFrame, filepath)
                residue_df, source_filepath = self._load_single_residue_file(row)
                
                if residue_df is None or not all(col in residue_df.columns for col in invariant_cols):
                    continue

                # MODIFICATION: Reset occurrence counter for each new file (row)
                occurrence_in_file_counter = 0

                for kmer in all_top_kmers:
                    if kmer not in row:
                        continue
                    
                    value = row[kmer]

                    if pd.api.types.is_scalar(value) and pd.isna(value):
                        continue

                    locations = []
                    try:
                        if isinstance(value, str) and value.startswith('['):
                            locations = ast.literal_eval(value)
                        elif isinstance(value, list):
                            locations = value
                    except (ValueError, SyntaxError):
                        continue

                    for start_loc in locations:
                        kmer_slice = residue_df.iloc[start_loc : start_loc + self.k]
                        if len(kmer_slice) == self.k:
                            for i, residue_row in kmer_slice.iterrows():
                                data_dict = {
                                    'position_in_kmer': i - start_loc,
                                    'source_file': source_filepath,
                                    'start_location': start_loc
                                }
                                for col in invariant_cols:
                                    data_dict[col] = residue_row[col]
                                temp_kmer_data_buffer[kmer].append(data_dict)
                            
                            # MODIFICATION: Increment counter after processing one full k-mer occurrence
                            occurrence_in_file_counter += 1
                
                processed_in_session.append(index)

                if len(processed_in_session) % checkpoint_interval == 0:
                    self._save_kmer_data(temp_kmer_data_buffer, output_dir)
                    with open(progress_file, 'a') as f_progress:
                        f_progress.write('\n'.join(map(str, processed_in_session)) + '\n')
                    
                    temp_kmer_data_buffer = {kmer: [] for kmer in all_top_kmers}
                    processed_in_session = []

        finally:
            print("\n--- Finalizing session: saving any remaining data... ---")
            self._save_kmer_data(temp_kmer_data_buffer, output_dir)
            if processed_in_session:
                with open(progress_file, 'a') as f_progress:
                    f_progress.write('\n'.join(map(str, processed_in_session)) + '\n')
            print("Save complete.")

        return self._load_final_data(output_dir)

    def _load_final_data(self, output_dir: str) -> Dict[str, pd.DataFrame]:
        """Helper to load all generated k-mer CSVs into a dictionary of DataFrames."""
        final_dataframes = {}
        all_top_kmers = [col for col in self.df_with_locations.columns if len(col) == self.k and isinstance(col, str)]
        for kmer in all_top_kmers:
            filepath = os.path.join(output_dir, f"{kmer}.csv")
            if os.path.exists(filepath):
                try:
                    final_dataframes[kmer] = pd.read_csv(filepath)
                except pd.errors.EmptyDataError:
                    final_dataframes[kmer] = pd.DataFrame() # Return empty df if file is empty
        return final_dataframes

import pandas as pd
import os
import glob
import ast
from collections import Counter
from typing import List, Dict, Tuple
from tqdm import tqdm

# This initializes tqdm for use with pandas .progress_apply()
tqdm.pandas(desc="Processing Rows")

# === Global Parameters ===
# Note: These Windows-style UNC paths will fail on Linux unless the remote share
# is mounted. For example, if \\rfs01\rdm01 is mounted at /mnt/rfs01, the path
# should be changed to /mnt/rfs01/GeometryOfProteins/...
GEO_DATA_PATH = "\\\\rfs01\\rdm01\\GeometryOfProteins\\Ziqiu\\invariant_result\\geo_info_18Mar2024.csv"
INVARIANT_DATA_PATH = "\\\\rfs01\\rdm01\\GeometryOfProteins\\Ziqiu\\invariant_result\\invariant_with_weak\\"

class KmerAnalyzer:
    """
    A class to handle the analysis of k-mers from protein geometric sequences.
    """
    def __init__(self, k: int, depth: int, geo_data_path: str = GEO_DATA_PATH, invariant_data_path: str = INVARIANT_DATA_PATH):
        """
        Initializes the analyzer, using global paths as defaults.
        """
        self.k = k
        self.depth = depth
        self.geo_data_path = geo_data_path
        self.invariant_data_path = invariant_data_path
        
        try:
            self.geo_df = pd.read_csv(self.geo_data_path)
            # Magic number -3 removes the trailing 'END' from sequences.
            self.sequences = self.geo_df['geo_seq'].dropna().astype(str).str[:-3].tolist()
        except FileNotFoundError:
            print(f"Error: The file '{self.geo_data_path}' was not found.")
            self.geo_df = pd.DataFrame()
            self.sequences = []
            
        self.df_with_locations = None
        self._file_cache = {}
        
        # --- OPTIMIZATION: Build the file map once during initialization ---
        self.invariant_file_map = {}
        self._create_invariant_file_map()

    def _create_invariant_file_map(self):
        """
        OPTIMIZATION: Scans the invariant data directory once and creates a mapping 
        from a unique identifier tuple to the full file path. This avoids repeated
        file system searches in the main processing loop.
        """
        print("Building invariant file map... this may take a moment.")
        search_path = os.path.join(self.invariant_data_path, "*.csv")
        all_files = glob.glob(search_path)

        for f_path in tqdm(all_files, desc="Mapping invariant files"):
            filename = os.path.basename(f_path)
            parts = filename.split('-')
            if len(parts) >= 4:
                try:
                    pdb_id = parts[0]
                    residue_id = parts[1]
                    model_id = parts[2]
                    chain_id = parts[3]
                    key = (pdb_id, residue_id, model_id, chain_id)
                    
                    if key not in self.invariant_file_map:
                        self.invariant_file_map[key] = f_path
                except IndexError:
                    pass
        print(f"File map built. Found {len(self.invariant_file_map)} unique invariant files.")

    def find_top_kmers(self) -> List[Tuple[str, int]]:
        """
        Finds the most frequent k-mers and their counts across all sequences.
        Returns a list of (k-mer, count) tuples.
        """
        all_kmers = []
        for s in tqdm(self.sequences, desc=f"Finding all {self.k}-mers"):
            if len(s) >= self.k:
                for i in range(len(s) - self.k + 1):
                    all_kmers.append(s[i:i+self.k])
        
        if not all_kmers:
            return []

        frequencies = Counter(all_kmers)
        top_kmers_with_counts = frequencies.most_common(self.depth)
        return top_kmers_with_counts

    def _find_single_kmer_locations(self, sequence: str, kmer_to_find: str) -> List[int]:
        """Helper function to find all locations of a specific k-mer in a sequence."""
        locations = []
        k = len(kmer_to_find)
        if not isinstance(sequence, str):
            return locations
            
        for i in range(len(sequence) - k + 1):
            if sequence[i:i+k] == kmer_to_find:
                locations.append(i)
        return locations

    def create_location_data(self, force_rerun: bool = False):
        """
        Creates and saves a DataFrame with location information for the top k-mers.
        """
        output_dir = f"k{self.k}"
        output_filename = os.path.join(output_dir, f"geo_info_k{self.k}_locations.csv")
        os.makedirs(output_dir, exist_ok=True)

        run_analysis = True
        if os.path.exists(output_filename) and not force_rerun:
            header_df = pd.read_csv(output_filename, nrows=0)
            kmer_columns_in_file = [col for col in header_df.columns if len(col) == self.k]
            
            if len(kmer_columns_in_file) >= self.depth:
                print(f"Loading existing data from '{output_filename}'...")
                self.df_with_locations = pd.read_csv(output_filename)
                run_analysis = False
            else:
                print(f"File only contains data for depth {len(kmer_columns_in_file)}. Re-running for depth {self.depth}...")

        if run_analysis:
            print("Running full analysis to generate location data...")
            top_kmers_with_counts = self.find_top_kmers()
            if not top_kmers_with_counts:
                self.df_with_locations = pd.DataFrame()
                return self.df_with_locations

            top_kmers = [kmer for kmer, count in top_kmers_with_counts]
            self.df_with_locations = self.geo_df.copy()
            
            for kmer in tqdm(top_kmers, desc="Finding k-mer locations"):
                self.df_with_locations[kmer] = self.df_with_locations['geo_seq'].dropna().astype(str).apply(
                    lambda seq: self._find_single_kmer_locations(sequence=seq, kmer_to_find=kmer)
                )
            
            mask = self.df_with_locations[top_kmers].apply(
                lambda row: any(isinstance(loc_list, list) and len(loc_list) > 0 for loc_list in row),
                axis=1
            )
            self.df_with_locations = self.df_with_locations[mask]

            self.df_with_locations.to_csv(output_filename, index=False)
            print(f"Saved new location data to '{output_filename}'.")
        
        return self.df_with_locations

    def _load_single_residue_file(self, row: pd.Series) -> Tuple[pd.DataFrame | None, str | None]:
        """
        MODIFIED: Finds and loads the invariant data file for a row using the 
        pre-built file map, with caching. 
        Returns a tuple of (DataFrame, filepath).
        """
        try:
            key = (
                row['pdb_id'],
                str(row['residue_id']),
                str(row['model_id']),
                row['chain_id']
            )
        except KeyError:
            return None, None

        filepath = self.invariant_file_map.get(key)

        if not filepath:
            return None, None

        if filepath in self._file_cache:
            # Return cached dataframe and the filepath
            return self._file_cache[filepath], filepath

        try:
            data_df = pd.read_csv(filepath)
            self._file_cache[filepath] = data_df 
            return data_df, filepath
        except Exception: 
            self._file_cache[filepath] = None 
            return None, None

    def _save_kmer_data(self, kmer_data_collection: Dict[str, List], output_dir: str):
        """Helper method to save the collected k-mer data to CSV files."""
        for kmer, data_list in kmer_data_collection.items():
            if data_list:
                filepath = os.path.join(output_dir, f"{kmer}.csv")
                write_header = not os.path.exists(filepath)
                df_to_save = pd.DataFrame(data_list)
                df_to_save.to_csv(filepath, mode='a', header=write_header, index=False)

    def extract_invariant_data(self, force_rerun: bool = False, checkpoint_interval: int = 100) -> Dict[str, pd.DataFrame]:
        """
        MODIFIED: Extracts the 9-part invariant data for each top k-mer.
        - occurrence_id is now local to each source file.
        - source_file column now contains the full path to the data file.
        """
        if self.df_with_locations is None:
            print("Location data not found. Running .create_location_data() first...")
            self.create_location_data()
            if self.df_with_locations is None or self.df_with_locations.empty:
                print("Failed to create location data. Aborting.")
                return {}

        output_dir = f"k{self.k}"
        os.makedirs(output_dir, exist_ok=True)
        progress_file = os.path.join(output_dir, '.progress.csv')
        
        df_to_process = self.df_with_locations.copy()
        processed_indices = set()

        if os.path.exists(progress_file) and not force_rerun:
            print(f"Found progress file: '{progress_file}'. Resuming analysis.")
            try:
                progress_df = pd.read_csv(progress_file, header=None, names=['processed_index']).dropna()
                processed_indices = set(progress_df['processed_index'].astype(int))
                df_to_process = df_to_process[~df_to_process.index.isin(processed_indices)]
            except (pd.errors.EmptyDataError, KeyError):
                print("Progress file is empty or malformed. Starting from scratch.")
                if os.path.exists(progress_file):
                    os.remove(progress_file)

        if df_to_process.empty:
            print("All rows have already been processed. Nothing to do.")
            # Still load and return existing data
            return self._load_final_data(output_dir)

        all_top_kmers = [col for col in self.df_with_locations.columns if len(col) == self.k and isinstance(col, str)]
        
        temp_kmer_data_buffer = {kmer: [] for kmer in all_top_kmers}
        
        invariant_cols = [
            'length(N)', 'length(A)', 'length(C)', 
            'angle(N)', 'angle(A)', 'angle(C)',
            'tau(NA)', 'tau(AC)', 'tau(CN)'
        ]

        processed_in_session = []
        
        try:
            for index, row in tqdm(df_to_process.iterrows(), total=len(df_to_process), desc="Processing files"):
                # MODIFICATION: Unpack tuple of (DataFrame, filepath)
                residue_df, source_filepath = self._load_single_residue_file(row)
                
                if residue_df is None or not all(col in residue_df.columns for col in invariant_cols):
                    continue

                # MODIFICATION: Reset occurrence counter for each new file (row)
                occurrence_in_file_counter = 0

                for kmer in all_top_kmers:
                    if kmer not in row:
                        continue
                    
                    value = row[kmer]

                    if pd.api.types.is_scalar(value) and pd.isna(value):
                        continue

                    locations = []
                    try:
                        if isinstance(value, str) and value.startswith('['):
                            locations = ast.literal_eval(value)
                        elif isinstance(value, list):
                            locations = value
                    except (ValueError, SyntaxError):
                        continue

                    for start_loc in locations:
                        kmer_slice = residue_df.iloc[start_loc : start_loc + self.k]
                        if len(kmer_slice) == self.k:
                            for i, residue_row in kmer_slice.iterrows():
                                data_dict = {
                                    'position_in_kmer': i - start_loc,
                                    'source_file': source_filepath,
                                    'start_location': start_loc
                                }
                                for col in invariant_cols:
                                    data_dict[col] = residue_row[col]
                                temp_kmer_data_buffer[kmer].append(data_dict)
                            
                            # MODIFICATION: Increment counter after processing one full k-mer occurrence
                            occurrence_in_file_counter += 1
                
                processed_in_session.append(index)

                if len(processed_in_session) % checkpoint_interval == 0:
                    self._save_kmer_data(temp_kmer_data_buffer, output_dir)
                    with open(progress_file, 'a') as f_progress:
                        f_progress.write('\n'.join(map(str, processed_in_session)) + '\n')
                    
                    temp_kmer_data_buffer = {kmer: [] for kmer in all_top_kmers}
                    processed_in_session = []

        finally:
            print("\n--- Finalizing session: saving any remaining data... ---")
            self._save_kmer_data(temp_kmer_data_buffer, output_dir)
            if processed_in_session:
                with open(progress_file, 'a') as f_progress:
                    f_progress.write('\n'.join(map(str, processed_in_session)) + '\n')
            print("Save complete.")

        return self._load_final_data(output_dir)

    def _load_final_data(self, output_dir: str) -> Dict[str, pd.DataFrame]:
        """Helper to load all generated k-mer CSVs into a dictionary of DataFrames."""
        final_dataframes = {}
        all_top_kmers = [col for col in self.df_with_locations.columns if len(col) == self.k and isinstance(col, str)]
        for kmer in all_top_kmers:
            filepath = os.path.join(output_dir, f"{kmer}.csv")
            if os.path.exists(filepath):
                try:
                    final_dataframes[kmer] = pd.read_csv(filepath)
                except pd.errors.EmptyDataError:
                    final_dataframes[kmer] = pd.DataFrame() # Return empty df if file is empty
        return final_dataframes
