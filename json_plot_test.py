import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Define the age range for snapshots 300 to 385
min_age = 6.58906279
max_age = 8.81250219

def create_custom_colormap():
    colors = [(0.6, 0.8, 1), (1, 0.6, 0.8)]  # Light Blue to Pink
    cmap_name = 'light_blue_to_pink'
    return mcolors.LinearSegmentedColormap.from_list(cmap_name, colors)

def age_formatter(x, pos):
    snapshot_index = int(x)
    age = min_age + (snapshot_index - 300) * (max_age - min_age) / (385 - 300)
    return f'{age:.2f}'

def plot_density_vs_radial_bins(json_file, save_path=None):
    # Read the data from the JSON file
    with open(json_file, 'r') as file:
        profiles = json.load(file)
    
    # Determine the range of snapshot indices
    min_snapshot = min(int(name.split('_')[1]) for name in profiles.keys())
    max_snapshot = max(int(name.split('_')[1]) for name in profiles.keys())
    num_snapshots = max_snapshot - min_snapshot + 1
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create the custom colormap
    custom_cmap = create_custom_colormap()
    
    # Initialize an array to store density values for each snapshot
    all_density_profiles = []
    
    # Loop through each snapshot in the profiles dictionary
    for snapshot_name, snapshot_data in profiles.items():
        pStar_data = snapshot_data['pStar']
        
        radial_bins = np.array(pStar_data['radial_bins_kpc'])
        density = np.array(pStar_data['density'])
        
        # Add the density profile to the array
        all_density_profiles.append(density)
        
        # Get the index for color from colormap based on snapshot number
        snapshot_index = int(snapshot_name.split('_')[1])
        
        # Normalize the snapshot_index to map to colormap range
        norm = Normalize(vmin=min_snapshot, vmax=max_snapshot)
        color = custom_cmap(norm(snapshot_index))
        
        ax.plot(radial_bins, density, color=color)
        #ax.semilogy()
        
    # Calculate and plot the median profile
    median_density_profile = np.median(all_density_profiles, axis=0)
    ax.plot(radial_bins, median_density_profile, color='black', linestyle='dashed', label='Median Profile')
    
    # Add labels and title
    ax.set_xlabel('R [kpc]')
    ax.set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
    ax.set_title('Stellar Density')
    
#     ax.plot(pStar['rbins'].in_units('kpc'),pStar['vr_disp'].in_units('km s^-1'))
#     ax.set_title('Stellar Radial Velocity Dispersion')
#     ax.set_xlabel('R [kpc]')
#     ax.set_ylabel('$\sigma_{r}$')
    
#     ax.plot(pStar['rbins'].in_units('kpc'),pStar['vt_disp'].in_units('km s^-1'))
#     ax.set_title('Stellar Tangential Velocity Dispersion')
#     ax.set_xlabel('R [kpc]')
#     ax.set_ylabel('$\sigma_{r}$')
    
#     ax.plot(pStar['rbins'],pStar['density'])
#     ax.semilogy()
#     ax.set_title('Stellar Density')
#     ax.set_xlabel('R [kpc]')
#     ax.set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
    
#     ax.plot(pDark['rbins'].in_units('kpc'),pDark['vr_disp'].in_units('km s^-1'))
#     ax.set_title('Dark Matter Radial Velocity Dispersion')
#     ax.set_xlabel('R [kpc]')
#     ax.set_ylabel('$\sigma_{r}$')
    
#     ax.plot(pDark['rbins'].in_units('kpc'),pDark['vt_disp'].in_units('km s^-1'))
#     ax.set_title('Dark Matter Tangential Velocity Dispersion')
#     ax.set_xlabel('R [kpc]')
#     ax.set_ylabel('$\sigma_{r}$')
    
#     ax.plot(pDark['rbins'],pDark['density'])
#     ax.semilogy()
#     ax.set_title('Dark Matter Density')
#     aax.set_xlabel('R [kpc]')
#     ax.set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
    
    
    
    # Add legend with specified location
    ax.legend(loc='upper right')
    
    # Create a ScalarMappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
    sm.set_array([])  # Need to set a dummy array for it to work
    
    # Position the colorbar horizontally under the graph
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.6)  # Adjust 'pad' as needed
    
    # Add colorbar using the ScalarMappable with age labels
    cbar = plt.colorbar(sm, cax=cax, orientation='horizontal', label='Age (Gyr)',
                        format=FuncFormatter(age_formatter))
    
    # Save the plot as an image or show it
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        print(f"Plot saved as '{save_path}'")
    else:
        plt.show()

# Call the function to save the plot as an image
plot_density_vs_radial_bins('profiles.json', save_path='density_plot.png')
