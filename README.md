### K-mer Invariant Data (`<kmer>.csv`)

For each top k-mer, a separate CSV file is created (e.g., `AAL.csv`). Each row in this file represents a single residue within an occurrence of that k-mer, detailing its geometric properties.

-----

\<img width="1858" height="817" alt="image" src="[[A_ramachandran.gif]]" /\>

**Columns:**

  * `position_in_kmer`: The 0-indexed position of the residue within the k-mer itself (from `0` to `k-1`).
  * `source_file`: The full filesystem path to the source `.csv` file from which the invariant data for this residue was extracted.
  * `start_location`: The starting index of this specific k-mer occurrence within the `source_file`.
  * `length(N)`, `length(A)`, `length(C)`: The three length invariants for the residue.
  * `angle(N)`, `angle(A)`, `angle(C)`: The three angle invariants for the residue.
  * `tau(NA)`, `tau(AC)`, `tau(CN)`: The three torsional angle (tau) invariants for the residue.

### Location Data (`seq_info_k{k}_locations.csv`)

This file acts as an intermediate lookup table. It is a filtered version of the input `seq_info...csv` file, containing only the rows (proteins) where at least one of the top k-mers was found. It adds new columns for each of the top k-mers.

## Class: `KmerAnalyzer`

This class contains all the logic for the analysis.

### Methods

  * `__init__(self, k: int, depth: int, ...)`
    The constructor initializes the analyzer. It sets the k-mer length **`k`** and the number of top k-mers to find **`depth`**. It loads the main protein sequence data and, for optimization, immediately scans the invariant data directory to build a file map. This map allows for rapid file lookups during the main processing step.

  * `find_top_kmers(self)`
    This method iterates through all protein sequences to count the occurrences of every possible k-mer. It returns a list of the **`depth`** most frequent k-mers and their counts.

  * `create_location_data(self, force_rerun: bool = False)`
    This method generates the `seq_info_k{k}_locations.csv` file described above. It uses the results from `find_top_kmers()` to find the specific locations of top k-mers in each protein sequence. To save time, it will load the data from an existing file if one is present, unless `force_rerun` is set to `True`.

  * `extract_invariant_data(self, force_rerun: bool = False, checkpoint_interval: int = 100)`
    This is the main data extraction method. It iterates through the location data file, and for each k-mer occurrence, it:

    1.  Loads the corresponding detailed invariant data file using the pre-built map.
    2.  Extracts the 9 geometric invariants for each of the `k` residues in the occurrence.
    3.  Saves the data to the appropriate `<kmer>.csv` file.

    This method supports **resumption**. It logs its progress and can pick up where it left off if interrupted. Data is saved periodically based on the `checkpoint_interval`.

-----

# K-mer Visualizer

A command-line tool to generate interactive 3D Ramachandran-like surface plots from pre-processed k-mer data files (e.g., from `k5/GSSGS.csv`). It visualizes the distribution of φ and ψ torsion angles for a given k-mer or a batch of k-mers.

### Dependencies

  * `pandas`
  * `plotly`
  * `tqdm`
  * Optional (for GIF generation): `imageio`, `kaleido`

### Usage

The script can be run on a single k-mer or in batch mode on all k-mers of a given length.

```bash
python kmer_visualizer.py <target> [options]
```

**`<target>`**

  * **A k-mer string** (e.g., `GSSGS`) to run in single mode.
  * **An integer** (e.g., `5`) to run in batch mode for all k-mers of that length (e.g., all files in the `k5/` directory).

**`[options]`**

  * `--use_360`: Use a 0-360° range for angles instead of the default -180° to 180°.
  * `--colorscale <name>`: Specify a Plotly colorscale. **Default is `Magma`**.
  * `--log_scale`: Apply a logarithmic scale to both the z-axis (height) and the color map. Ideal for data with large frequency outliers.
  * `--gif`: Generate a rotating GIF instead of an interactive plot.
  * `--gif_name <filename>`: (Single k-mer mode only) Specify a custom output filename for the GIF.

### Examples

```bash
# Generate a standard interactive plot for 'GSSGS'
python kmer_visualizer.py GSSGS

# Generate a plot with a log scale and the 'Plasma' colorscale
python kmer_visualizer.py GSSGS --log_scale --colorscale Plasma

# Generate a rotating GIF file named 'GSSGS_ramachandran.gif'
python kmer_visualizer.py GSSGS --gif

# Run in batch mode to generate GIFs for all 5-mers found in the k5/ directory
python kmer_visualizer.py 5 --gif
```