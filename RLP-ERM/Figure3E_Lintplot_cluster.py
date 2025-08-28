import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set a modern style for matplotlib
plt.style.use('seaborn-v0_8-white')

# Load the merged data file
merged_file = "/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/Figures/ERM_gene_heatmap_cluster.tsv"
merged_data = pd.read_csv(merged_file, sep="\t", index_col=0)

# Define time points and their numeric values
time_points = ['6h', '12h', '24h', '48h', '72h']
time_numeric = [6, 12, 24, 48, 72]

# Get unique cluster labels (sorted)
clusters = sorted(merged_data['cluster'].unique())

# Create a cluster-to-color mapping
cluster_color_map = {
        1: "#E27508",   # brown
        2: "#F72FA0",  # red
        3: "#1996F6",  # blue
        4: "#963CFD",  # purple 
        5: "#FCB715",  # orange
        6: "#2CDC5B",  # green
        7: "#DED225"  # yellow
    }

# Loop through each cluster and plot
for cl in clusters:
    cluster_data = merged_data[merged_data['cluster'] == cl]
    
    # Create figure with smaller size
    fig, ax = plt.subplots(figsize=(2.5, 1.8))
    
    # Plot individual gene profiles in light gray with low alpha
    for gene, row in cluster_data.iterrows():
        ax.plot(time_numeric, row[time_points], color='#D3D3D3', alpha=0.1, linewidth=0.4)
    
    # Calculate and plot the mean profile with custom color
    cluster_mean = cluster_data[time_points].mean()
    ax.plot(time_numeric, cluster_mean, marker='o', color=cluster_color_map.get(cl, '#000000'), 
            linewidth=1.2, markersize=3, label=f"Cluster {cl} Mean")  # Default to black if cluster not in map
    
    # Customize axes
    ax.set_xlabel("Time (h)", fontsize=7, fontfamily='Arial')
    ax.set_ylabel("Editing", fontsize=7, fontfamily='Arial')
    ax.set_title(f"Cluster {cl}", fontsize=8, fontfamily='Arial', fontweight='bold', pad=4)
    
    # Set ticks explicitly
    ax.set_xticks(time_numeric)
    ax.set_yticks([-1, 0, 1])  # Only show -1, 0, 1 on y-axis
    
    # Adjust ticks to point outside with explicit length
    ax.tick_params(axis='x', which='major', direction='out', length=3, labelsize=5, pad=1)
    ax.tick_params(axis='y', which='major', direction='out', length=3, labelsize=5, pad=1)
    ax.tick_params(axis='y', which='minor', left=False)  # Hide minor y-axis ticks
    
    # Set frame (spines) to gray
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(0.4)
    
    # Set a light background
    fig.patch.set_facecolor('#F8F8F8')
    ax.set_facecolor('#FFFFFF')
    
    # Adjust layout to prevent clipping
    plt.tight_layout()
    
    # Save and show the plot
    plt.savefig(f"/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/Figures/cluster_{cl}_trend.png", 
                dpi=600, bbox_inches='tight')
    plt.show()