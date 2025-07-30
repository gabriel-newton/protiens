import numpy as np
import plotly.graph_objects as go

# --- Configuration ---
# A lower threshold is better, as the log scale will handle the differences.
COUNT_THRESHOLD = 20

def plot_3d_surface_log_color(input_npy="angle_density_grid.npy"):
    """
    Plots a 3D surface with a logarithmic z-axis AND a logarithmic
    color scale to visualize data with extreme outliers.
    """
    print(f"Loading grid from '{input_npy}'...")
    try:
        density_grid = np.load(input_npy)
    except FileNotFoundError:
        print(f"Error: The file '{input_npy}' was not found.")
        return

    # Filter data below the threshold by replacing with NaN
    surface_data = density_grid.astype(float)
    surface_data[surface_data <= COUNT_THRESHOLD] = np.nan
    
    # --- THIS IS THE FIX ---
    # Create a separate array for the color values by taking the log10 of the counts.
    # We add a very small number to avoid log(0) errors if any zeros remain.
    log_color_data = np.log10(surface_data + 1e-6)
    
    # Define the tick values and labels for the color bar
    cbar_ticks = [np.log10(100), np.log10(1000), np.log10(10000), np.log10(100000)]
    cbar_tick_text = ['100', '1k', '10k', '100k']
    # -----------------------

    angle_range = np.arange(-180, 180)
    
    print("Generating 3D surface plot...")
    
    fig = go.Figure(data=[
        go.Surface(
            z=surface_data,          # Z-axis height uses original data
            x=angle_range, 
            y=angle_range,
            surfacecolor=log_color_data, # Color is determined by the log-transformed data
            colorscale='jet',
            colorbar=dict(
                title='Count',
                tickvals=cbar_ticks,    # Set custom tick values
                ticktext=cbar_tick_text # Set custom tick labels
            )
        )
    ])
    
    fig.update_layout(
        title='3D Surface Plot of 2-mer Angle Density (Log Height & Log Color)',
        scene=dict(
            xaxis_title='Average Phi Angle (°)',
            yaxis_title='Average Psi Angle (°)',
            zaxis_title='Count (log scale)',
            zaxis_type='log' # Keep the Z-axis height logarithmic
        ),
        margin=dict(l=40, r=40, b=40, t=80)
    )

    print("Opening plot in browser...")
    fig.show(renderer="browser")

if __name__ == "__main__":
    plot_3d_surface_log_color()