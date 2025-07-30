<<<<<<< HEAD
import numpy as np
import pandas as pd
import os
import imageio
import shutil
import plotly.graph_objects as go
from tqdm import tqdm
import sys
import argparse

class KmerVisualizer:
    """
    A class to load k-mer data and generate interactive 3D visualizations
    of Ramachandran plots using Plotly.
    """
    def __init__(self, kmer: str):
        self.kmer = kmer
        k_value = len(kmer)
        self.input_file = os.path.join(f"k{k_value}", f"{kmer}.csv")
        self.df = None
        
        try:
            self.df = pd.read_csv(self.input_file)
            if __name__ != "__main__" or not sys.argv[1].isdigit():
                 print(f"Successfully loaded data for k-mer '{self.kmer}' from '{self.input_file}'.")
        except FileNotFoundError:
            raise FileNotFoundError(f"Error: The file '{self.input_file}' was not found.")
        except pd.errors.EmptyDataError:
            print(f"Warning: The file '{self.input_file}' is empty.")
            self.df = pd.DataFrame()

    def _build_figure(self, use_360_range: bool = False, colorscale: str = 'Magma', use_log_scale: bool = False) -> go.Figure:
        if self.df is None or self.df.empty:
            return None

        required_cols = ['tau(NA)', 'tau(AC)']
        if not all(col in self.df.columns for col in required_cols):
            return None

        angle_data = self.df[required_cols].dropna().copy()
        
        if use_360_range:
            angle_data['phi'] = (angle_data['tau(NA)'].round() % 360).astype(int)
            angle_data['psi'] = (angle_data['tau(AC)'].round() % 360).astype(int)
            full_range = np.arange(0, 361)
            axis_config = dict(range=[0, 360], tickmode='linear', dtick=60)
        else:
            angle_data['phi'] = angle_data['tau(NA)'].round().astype(int)
            angle_data['psi'] = angle_data['tau(AC)'].round().astype(int)
            full_range = np.arange(-180, 181)
            axis_config = dict(range=[-180, 180], tickmode='linear', dtick=60)
            
        freq_counts = angle_data.groupby(['phi', 'psi']).size().reset_index(name='count')
        z_data = freq_counts.pivot_table(index='psi', columns='phi', values='count', fill_value=0)
        z_data = z_data.reindex(index=full_range, columns=full_range, fill_value=0)
        
        # --- MODIFIED: Prepare z_data for both linear and log scales ---
        # Convert to float to allow for NaN values
        z_data_numeric = z_data.values.astype(float)
        z_data_numeric[z_data_numeric == 0] = np.nan
        # ----------------------------------------------------------------

        fig = go.Figure(data=[go.Surface(
            z=z_data_numeric, 
            x=z_data.columns, 
            y=z_data.index,
            colorscale=colorscale, 
            cmin=1,
            colorbar=dict(title='Residue Count', len=0.75), 
            connectgaps=False
        )])
        
        fig.update_layout(
            title=f"3D Ramachandran Plot for k-mer: '{self.kmer}'",
            scene=dict(
                xaxis_title='Phi (φ) / tau(NA) [degrees]',
                yaxis_title='Psi (ψ) / tau(AC) [degrees]',
                zaxis_title='Frequency Count',
                xaxis=axis_config, yaxis=axis_config,
            ),
            width=800, height=800, margin=dict(l=65, r=50, b=65, t=90)
        )
        
        if use_log_scale:
            # --- NEW: Apply log scale to both height and color map ---
            with np.errstate(invalid='ignore'): # Ignore warnings for log10(nan)
                log_color_data = np.log10(z_data_numeric)

            # Dynamically create ticks for the color bar (e.g., 10, 100, 1k)
            min_log = np.nanmin(log_color_data)
            max_log = np.nanmax(log_color_data)
            tick_vals = np.arange(np.ceil(min_log), np.floor(max_log) + 1)
            tick_text = [f"{10**v:.0f}" for v in tick_vals]
            
            # Update the surface trace to use the log-transformed color data
            fig.update_traces(
                surfacecolor=log_color_data,
                colorbar=dict(
                    title='Residue Count',
                    len=0.75,
                    tickvals=tick_vals,
                    ticktext=tick_text
                )
            )
            
            # Update the scene for a log z-axis (height)
            fig.update_scenes(
                zaxis_type='log',
                zaxis_title_text='Frequency Count (log scale)'
            )
            # -------------------------------------------------------------
            
        return fig
            
    def plot_3d_ramachandran(self, use_360_range: bool = False, colorscale: str = 'Magma', use_log_scale: bool = False):
        fig = self._build_figure(use_360_range=use_360_range, colorscale=colorscale, use_log_scale=use_log_scale)
        if fig:
            fig.show(renderer="browser")
    
    def create_ramachandran_gif(self, output_filename: str = "ramachandran_3d.gif", duration: float = 0.1, steps: int = 50, use_360_range: bool = False, colorscale: str = 'Magma', use_log_scale: bool = False):
        fig = self._build_figure(use_360_range=use_360_range, colorscale=colorscale, use_log_scale=use_log_scale)
        if not fig:
            return

        temp_dir = f"temp_gif_frames_{self.kmer}"
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
            
        print(f"Generating {steps} frames for '{self.kmer}' GIF...")
        frames = []
        for i in tqdm(range(steps), desc=f"Frames for {self.kmer}", leave=False):
            angle = (i / steps) * 2 * np.pi
            camera_eye = dict(x=2 * np.cos(angle), y=2 * np.sin(angle), z=1.5)
            fig.update_layout(scene_camera=dict(eye=camera_eye))
            frame_path = os.path.join(temp_dir, f"frame_{i:03d}.png")
            fig.write_image(frame_path)
            frames.append(imageio.imread(frame_path))
            
        print(f"Compiling frames into '{output_filename}'...")
        imageio.mimsave(output_filename, frames, duration=duration)
        shutil.rmtree(temp_dir)
        print(f"GIF creation for '{self.kmer}' complete.")

def process_kmer(kmer_name, args):
    try:
        visualizer = KmerVisualizer(kmer=kmer_name)
        if visualizer.df is not None and not visualizer.df.empty:
            if args.gif:
                output_filename = args.gif_name if args.gif_name and not args.target.isdigit() else f"{kmer_name}_ramachandran.gif"
                visualizer.create_ramachandran_gif(
                    output_filename=output_filename,
                    use_360_range=args.use_360, colorscale=args.colorscale, use_log_scale=args.log_scale
                )
            else:
                visualizer.plot_3d_ramachandran(
                    use_360_range=args.use_360, colorscale=args.colorscale, use_log_scale=args.log_scale
                )
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while processing '{kmer_name}': {e}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate 3D Ramachandran plots for k-mers from CSV data.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("target", type=str, help="A k-mer string or a k-value for batch processing.")
    parser.add_argument("--use_360", action="store_true", help="Use 0-360 degree range instead of -180 to 180.")
    parser.add_argument("--colorscale", type=str, default="Magma", help="Plotly colorscale to use.")
    parser.add_argument("--gif", action="store_true", help="Generate a rotating GIF instead of an interactive plot.")
    parser.add_argument("--gif_name", type=str, default=None, help="Optional: specify output filename for the GIF.")
    parser.add_argument("--log_scale", action="store_true", help="Apply a logarithmic scale to the frequency (z) axis and color map.")
    args = parser.parse_args()

    if args.target.isdigit():
        k_value = int(args.target)
        kmer_dir = f"k{k_value}"
        if not os.path.isdir(kmer_dir):
            print(f"Error: Directory '{kmer_dir}' not found.", file=sys.stderr)
            sys.exit(1)
        csv_files = sorted([f for f in os.listdir(kmer_dir) if f.lower().endswith('.csv')])
        if not csv_files:
            print(f"No .csv files found in '{kmer_dir}'.")
            sys.exit(0)
        print(f"--- Batch Mode: Found {len(csv_files)} k-mers to process in '{kmer_dir}' ---")
        for filename in tqdm(csv_files, desc="Processing k-mers"):
            kmer_name = os.path.splitext(filename)[0]
            process_kmer(kmer_name.upper(), args)
    else:
        kmer = args.target.upper()
        if args.gif:
             print(f"--- Single Mode: Generating GIF for '{kmer}' ---")
        else:
             print(f"--- Single Mode: Generating plot for '{kmer}' ---")
        process_kmer(kmer, args)
=======
import numpy as np
import pandas as pd
import os
import imageio
import shutil
import plotly.graph_objects as go
from tqdm import tqdm
import sys

class KmerVisualizer:
    """
    A class to load k-mer data and generate interactive 3D visualizations
    of Ramachandran plots using Plotly.
    """
    def __init__(self, kmer: str):
        """
        Initializes the visualizer by loading the data for a specific k-mer.

        Args:
            kmer (str): The k-mer string to analyze (e.g., 'AAAAA').
        """
        self.kmer = kmer
        k_value = len(kmer)
        self.input_file = os.path.join(f"k{k_value}", f"{kmer}.csv")
        self.df = None
        
        try:
            self.df = pd.read_csv(self.input_file)
            print(f"Successfully loaded data for k-mer '{self.kmer}' from '{self.input_file}'.")
        except FileNotFoundError:
            print(f"Error: The file '{self.input_file}' was not found.")
            sys.exit(1) # Exit if the file is not found
        except pd.errors.EmptyDataError:
            print(f"Warning: The file '{self.input_file}' is empty.")

    def _build_figure(self, use_360_range: bool = False, colorscale: str = 'Viridis') -> go.Figure:
        """
        Private helper method to prepare data and construct the Figure object.
        
        Args:
            use_360_range (bool): If True, display angles from 0-360 degrees.
                                  Otherwise, use -180 to 180 degrees.
            colorscale (str): The name of the Plotly colorscale to use.
        """
        if self.df is None or self.df.empty:
            print("Cannot generate plot: DataFrame is not loaded or is empty.")
            return None

        required_cols = ['tau(NA)', 'tau(AC)']
        if not all(col in self.df.columns for col in required_cols):
            print(f"Error: DataFrame must contain {required_cols} columns.")
            return None

        angle_data = self.df[required_cols].dropna().copy()
        
        # Toggle angle range based on the parameter
        if use_360_range:
            angle_data['phi'] = (angle_data['tau(NA)'].round() % 360).astype(int)
            angle_data['psi'] = (angle_data['tau(AC)'].round() % 360).astype(int)
            full_range = np.arange(0, 361)
            axis_config = dict(range=[0, 360], tickmode='linear', dtick=60)
        else:
            angle_data['phi'] = angle_data['tau(NA)'].round().astype(int)
            angle_data['psi'] = angle_data['tau(AC)'].round().astype(int)
            full_range = np.arange(-180, 181)
            axis_config = dict(range=[-180, 180], tickmode='linear', dtick=60)
            
        freq_counts = angle_data.groupby(['phi', 'psi']).size().reset_index(name='count')
        
        z_data = freq_counts.pivot_table(index='psi', columns='phi', values='count', fill_value=0)
        
        # Reindex to ensure a full grid for the selected range
        z_data = z_data.reindex(index=full_range, columns=full_range, fill_value=0)

        # Replace 0s with NaN to make the floor transparent.
        z_data[z_data == 0] = np.nan

        fig = go.Figure(data=[go.Surface(
            z=z_data.values,
            x=z_data.columns,
            y=z_data.index,
            colorscale=colorscale,
            cmin=1,
            colorbar=dict(title='Residue Count', len=0.75),
            connectgaps=False # Ensure gaps are not filled
        )])
        
        fig.update_layout(
            title=f"3D Ramachandran Plot for k-mer: '{self.kmer}'",
            scene=dict(
                xaxis_title='Phi (φ) / tau(NA) [degrees]',
                yaxis_title='Psi (ψ) / tau(AC) [degrees]',
                zaxis_title='Frequency Count',
                xaxis=axis_config,
                yaxis=axis_config,
            ),
            width=800,
            height=800,
            margin=dict(l=65, r=50, b=65, t=90)
        )
        return fig
            
    def plot_3d_ramachandran(self, use_360_range: bool = False, colorscale: str = 'Viridis'):
        """
        Creates and displays an interactive 3D Ramachandran surface plot.

        Args:
            use_360_range (bool): If True, display angles from 0-360 degrees.
                                  Defaults to False (-180 to 180 degrees).
            colorscale (str): The name of the Plotly colorscale to use.
        """
        fig = self._build_figure(use_360_range=use_360_range, colorscale=colorscale)
        if fig:
            # Explicitly set the renderer to open in a web browser.
            fig.show(renderer="browser")
    
    def create_ramachandran_gif(self, output_filename: str = "ramachandran_3d.gif", duration: float = 0.1, steps: int = 50, use_360_range: bool = False, colorscale: str = 'Viridis'):
        """
        Creates a rotating GIF of the 3D Ramachandran plot.

        Args:
            output_filename (str): The name of the output GIF file.
            duration (float): The duration (in seconds) of each frame in the GIF.
            steps (int): The number of frames to generate for the rotation.
            use_360_range (bool): If True, display angles from 0-360 degrees.
                                  Defaults to False (-180 to 180 degrees).
            colorscale (str): The name of the Plotly colorscale to use.
        """
        fig = self._build_figure(use_360_range=use_360_range, colorscale=colorscale)
        if not fig:
            return

        temp_dir = "temp_gif_frames"
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)
            
        print(f"Generating {steps} frames for GIF...")
        
        frames = []
        for i in tqdm(range(steps), desc="Creating GIF frames"):
            angle = (i / steps) * 2 * np.pi
            camera_eye = dict(x=2 * np.cos(angle), y=2 * np.sin(angle), z=1.5)
            
            fig.update_layout(scene_camera=dict(eye=camera_eye))
            
            frame_path = os.path.join(temp_dir, f"frame_{i:03d}.png")
            fig.write_image(frame_path)
            frames.append(imageio.imread(frame_path))
            
        print(f"Compiling frames into '{output_filename}'...")
        imageio.mimsave(output_filename, frames, duration=duration)
        
        # Clean up temporary files
        shutil.rmtree(temp_dir)
        print("GIF creation complete.")

if __name__ == "__main__":
    # Check for the correct number of arguments
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        print("Usage: python kmer_visualizer.py <KMER_STRING> [origin] [COLORMAP]")
        print("\nExamples:")
        print("  python kmer_visualizer.py AAAAA")
        print("  python kmer_visualizer.py GGGGG origin")
        print("  python kmer_visualizer.py CCCCC Plasma")
        print("  python kmer_visualizer.py TTTTT origin Cividis")
        sys.exit(1)

    kmer_arg = sys.argv[1].upper()
    
    # Set defaults for optional arguments
    use_360 = False
    colormap = 'Viridis'

    # Process optional arguments based on their position
    if len(sys.argv) >= 3:
        # Check if the third argument is 'origin'
        if sys.argv[2].lower() == 'origin':
            use_360 = True
            # If so, the fourth argument (if it exists) is the colormap
            if len(sys.argv) == 4:
                colormap = sys.argv[3]
        else:
            # If the third argument is not 'origin', assume it's the colormap
            colormap = sys.argv[2]
            # Ensure the fourth argument is not present in this case
            if len(sys.argv) == 4:
                print("Error: Cannot specify a colormap as both the 3rd and 4th argument.")
                sys.exit(1)

    # Instantiate the visualizer with the k-mer from the command line
    visualizer = KmerVisualizer(kmer=kmer_arg)
    
    # Check if data was loaded successfully before attempting to plot
    if visualizer.df is not None and not visualizer.df.empty:
        # Generate and display the interactive 3D plot with specified options
        print(f"Generating plot with 360-range={use_360} and colormap='{colormap}'")
        visualizer.plot_3d_ramachandran(use_360_range=use_360, colorscale=colormap)
>>>>>>> c528bcc31578f958837be719ab09982aba2642b4
