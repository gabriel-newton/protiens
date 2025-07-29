# K-mer Invariant Analyzer

This script analyzes geometric sequences of proteins to identify frequently occurring k-mers (subsequences of length *k*) and extracts their corresponding geometric invariant data. It processes a main data file of protein sequences, finds the most common k-mers, locates their occurrences, and then pulls detailed invariant measurements for each residue in those occurrences from a larger dataset of individual residue files.

---

### K-mer Invariant Data (`<kmer>.csv`)

For each top k-mer, a separate CSV file is created (e.g., `AAL.csv`). Each row in this file represents a single residue within an occurrence of that k-mer, detailing its geometric properties.

**Columns:**

* `position_in_kmer`: The 0-indexed position of the residue within the k-mer itself (from `0` to `k-1`).
* `source_file`: The full filesystem path to the source `.csv` file from which the invariant data for this residue was extracted.
* `start_location`: The starting index of this specific k-mer occurrence within the `source_file`.
* `length(N)`, `length(A)`, `length(C)`: The three length invariants for the residue.
* `angle(N)`, `angle(A)`, `angle(C)`: The three angle invariants for the residue.
* `tau(NA)`, `tau(AC)`, `tau(CN)`: The three torsional angle (tau) invariants for the residue.


### Location Data (`geo_info_k{k}_locations.csv`)

This file acts as an intermediate lookup table. It is a filtered version of the input `geo_info...csv` file, containing only the rows (proteins) where at least one of the top k-mers was found. It adds new columns for each of the top k-mers.


## Class: `KmerAnalyzer`

This class contains all the logic for the analysis.

### Methods

* `__init__(self, k: int, depth: int, ...)`
    The constructor initializes the analyzer. It sets the k-mer length **`k`** and the number of top k-mers to find **`depth`**. It loads the main protein sequence data and, for optimization, immediately scans the invariant data directory to build a file map. This map allows for rapid file lookups during the main processing step.

* `find_top_kmers(self)`
    This method iterates through all geometric sequences to count the occurrences of every possible k-mer. It returns a list of the **`depth`** most frequent k-mers and their counts.

* `create_location_data(self, force_rerun: bool = False)`
    This method generates the `geo_info_k{k}_locations.csv` file described above. It uses the results from `find_top_kmers()` to find the specific locations of top k-mers in each protein sequence. To save time, it will load the data from an existing file if one is present, unless `force_rerun` is set to `True`.

* `extract_invariant_data(self, force_rerun: bool = False, checkpoint_interval: int = 100)`
    This is the main data extraction method. It iterates through the location data file, and for each k-mer occurrence, it:
    1.  Loads the corresponding detailed invariant data file using the pre-built map.
    2.  Extracts the 9 geometric invariants for each of the `k` residues in the occurrence.
    3.  Saves the data to the appropriate `<kmer>.csv` file.

    This method supports **resumption**. It logs its progress and can pick up where it left off if interrupted. Data is saved periodically based on the `checkpoint_interval`.

# K-mer Visualizer

This script generates an interactive 3D Ramachandran plot for a given k-mer. It expects input data from a file path like k5/AAAAA.csv.

Dependencies

```pandas, plotly, tqdm, imageio, kaleido```

## Usage

python kmer_visualizer.py <KMER> [origin] [COLORMAP]

> <KMER>: (Required) The k-mer sequence (e.g., AAAAA).
>
> [origin]: (Optional) Use origin for a 0-360° plot range.
>
> [COLORMAP]: (Optional) A Plotly colorscale name (e.g., Plasma).

Examples
```
# Default plot
python kmer_visualizer.py AAAAA

# Plot with 360° range and a custom colormap
python kmer_visualizer.py GGGGG origin Plasma
```
