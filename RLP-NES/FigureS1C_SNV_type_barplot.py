import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_paths = [
    '/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_APO/APO_NES_SNV_summary.tsv',
    '/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_APO/APO_NES_3PABP_SNV_summary.tsv',
    '/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/TadA_NES_SNV_summary.tsv', 
    '/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/TadA_NES_3PABP_SNV_summary.tsv'
]

# Labels for the legend
file_labels = ['APO-NES', 'RLP-NES(APO)', 'TadA-NES', 'RLP-NES(TadA)'] 

# Colors for each file
colors = ['#9A9797', '#A9B31F', '#9F801A', '#E4317F']

# Read and process data from all files
data_frames = []
for i, file_path in enumerate(file_paths):
    try:
        # Read the file
        df = pd.read_csv(file_path, sep='\t')
        
        # Calculate percentage of each mutation type
        total_count = df['editing_ratio_count'].sum()
        df['percentage'] = (df['editing_ratio_count'] / total_count) * 100
        
        # Add file identifier
        df['file_id'] = i
        df['file_label'] = file_labels[i]
        
        data_frames.append(df)
        print(f"Successfully loaded {file_path}: {len(df)} mutation types, total count: {total_count}")
        
    except FileNotFoundError:
        print(f"Warning: File {file_path} not found. Please update the file path.")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

if not data_frames:
    print("No files were successfully loaded. Please check file paths.")
    exit()

# Combine all data
combined_df = pd.concat(data_frames, ignore_index=True)

# Get all unique mutation types and sort them
mutation_types = sorted(combined_df['mutation_type'].unique())
print(f"Found mutation types: {mutation_types}")

# Create the plot
fig, ax = plt.subplots(figsize=(4, 3.0), dpi=300)

# Set up bar positions
n_files = len(data_frames)
n_mutations = len(mutation_types)
bar_width = 0.8 / n_files
x_positions = np.arange(n_mutations)

# Create bars for each file
for i, (df, label, color) in enumerate(zip(data_frames, file_labels, colors)):
    # Create a series with all mutation types, filling missing ones with 0
    mutation_percentages = []
    for mut_type in mutation_types:
        matching_rows = df[df['mutation_type'] == mut_type]
        if len(matching_rows) > 0:
            mutation_percentages.append(matching_rows['percentage'].iloc[0])
        else:
            mutation_percentages.append(0)
    
    # Calculate bar positions for this file
    bar_positions = x_positions + (i - (n_files-1)/2) * bar_width
    
    # Create bars
    bars = ax.bar(bar_positions, mutation_percentages, bar_width, 
                  label=label, color=color, alpha=0.7)


# Customize the plot
ax.set_xlabel('Mutation Type', fontsize=10, fontweight='bold')
ax.set_ylabel('Percentage of Total Mutation Count (%)', fontsize=10, fontweight='bold')

# Set x-axis labels
ax.set_xticks(x_positions)
ax.set_xticklabels(mutation_types, rotation=45, ha='right')

# Add legend
ax.legend(bbox_to_anchor=(0.5, 1), loc='upper left', frameon=False)

# Set y-axis to start from 0
ax.set_ylim(0, None)

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Display statistics
print("\nSummary Statistics:")
print("="*50)
for i, (df, label) in enumerate(zip(data_frames, file_labels)):
    total_count = df['editing_ratio_count'].sum()
    most_common = df.loc[df['editing_ratio_count'].idxmax()]
    print(f"{label}:")
    print(f"  Total mutations: {total_count:,}")
    print(f"  Most common: {most_common['mutation_type']} ({most_common['percentage']:.1f}%)")
    print()

plt.savefig('/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/Figures/mutation_types_barplot.png', dpi=600, bbox_inches='tight')
# print("Plot saved as 'mutation_types_barplot.png'")

# Show the plot
plt.show()