import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

def plot_density_vs_radial_bins(json_file):
    # Read the data from the JSON file
    with open(json_file, 'r') as file:
        profiles = json.load(file)
    
    # Get a colormap for the color gradient
    cmap = get_cmap('viridis')
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Loop through each snapshot in the profiles dictionary
    for snapshot_name, snapshot_data in profiles.items():
        pStar_data = snapshot_data['pStar']
        
        radial_bins = np.array(pStar_data['radial_bins_kpc'])
        density = np.array(pStar_data['density'])
        
        # Get the index for color from colormap based on snapshot number
        snapshot_index = int(snapshot_name.split('_')[1])
        color = cmap(snapshot_index / len(profiles))
        
        ax.plot(radial_bins, density, color=color, label=f'Snapshot {snapshot_index}')
    
    # Add labels and title
    ax.set_xlabel('Radial Bins (kpc)')
    ax.set_ylabel('Density')
    ax.set_title('Density vs Radial Bins for Different Snapshots')
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=len(profiles)))
    sm.set_array([])
    cbar = plt.colorbar(sm, label='Snapshot Number')
    
    # Add legend
    ax.legend()
    
    # Show the plot
    plt.show()
    
# Call the function to plot the graph
plot_density_vs_radial_bins('profiles.json')
