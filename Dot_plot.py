import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.cluster import DBSCAN

# Configuration
base_path = "<path_to_input_dir>"
output_path = "<path_to_output_dir>"

populations = ["PJL", "BEB", "ITU"]
coverages = ["6x", "9x", "12x", "15x"]
colors = {"PJL": "darkred", "BEB": "darkgreen", "ITU": "navy"}

marker_size_small = 6
marker_size_12x = 12
linewidth_12x = 1.2

cluster_eps = 2e7
min_samples = 3

# Slot index for Y-axis offset
def slot_index(pop, cov):
    if cov == "6x":
        mapping = {"PJL": 3, "BEB": 2, "ITU": 1.0}
    elif cov == "9x":
        mapping = {"PJL": 3, "BEB": 2, "ITU": 1}
    elif cov == "12x":
        mapping = {"PJL": 3, "BEB": 2, "ITU": 1}
    else:  # 15x
        mapping = {"PJL": 3, "BEB": 2, "ITU": 1}
    return mapping[pop]

# Detect chromosome arms using DBSCAN
def detect_chromosome_arms(df, eps=cluster_eps, min_samples=min_samples):
    if len(df) < min_samples:
        return df
    
    try:
        clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(df[['X']])
        df['cluster'] = clustering.labels_
        df_clean = df[df['cluster'] != -1].copy()
        
        if len(df_clean) == 0:
            return df
        
        if df_clean['cluster'].nunique() > 1:
            cluster_info = []
            for cluster_id in df_clean['cluster'].unique():
                cluster_data = df_clean[df_clean['cluster'] == cluster_id]
                median_x = cluster_data['X'].median()
                cluster_info.append((cluster_id, median_x))
            
            cluster_info.sort(key=lambda x: x[1])
            
            merged_dfs = []
            for cluster_id, _ in cluster_info:
                cluster_df = df_clean[df_clean['cluster'] == cluster_id].copy()
                cluster_df = cluster_df.sort_values('X')
                merged_dfs.append(cluster_df[['X', 'Y']])
            
            result_df = pd.concat(merged_dfs, ignore_index=True)
        else:
            result_df = df_clean[['X', 'Y']].copy()
        
        return result_df.sort_values('X').reset_index(drop=True)
    
    except Exception as e:
        print(f"Warning: Arm detection failed, using original data: {e}")
        return df[['X', 'Y']].copy()

# Fix orientation and slope
def fix_orientation_and_slope(df):
    if len(df) < 2:
        return df
    
    # Ensure left-to-right orientation
    if df["X"].iloc[0] > df["X"].iloc[-1]:
        df = df.iloc[::-1].reset_index(drop=True)
    
    df = df.sort_values("X").reset_index(drop=True)
    
    # Flip Y-axis if slope is negative
    if len(df) >= 2:
        slope = np.polyfit(df["X"], df["Y"], 1)[0]
        if slope < 0:
            y_mid = (df["Y"].min() + df["Y"].max()) / 2
            df["Y"] = 2 * y_mid - df["Y"]
    
    return df
 
# Apply Y-axis offset for visualization (modify offset_factor as per chromosome)
def offset_data(df, pop, cov, y_min, y_max):
    df = df.copy()
    slot = slot_index(pop, cov)
    y_range = max(y_max - y_min, 1e-9)
    
    if cov in ["6x", "9x"]:
        offset_factor = 1.5e7
        df["Y_norm"] = (df["Y"] - y_min) / y_range
        df["Y_offset"] = df["Y_norm"] * 0.6e8 + slot * offset_factor
    elif cov == "12x":
        offset_factor = 0.3e8
        df["Y_norm"] = (df["Y"] - y_min) / y_range
        df["Y_offset"] = df["Y_norm"] * 0.6e8 + slot * offset_factor
    else:  # 15x
        offset_factor = 0.7e8
        df["Y_offset"] = df["Y"] + slot * offset_factor
    
    return df

# Main plotting function
def plot_chromosome_quadraplot(chr_name):
    os.makedirs(output_path, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True, dpi=300)
    axes = axes.flatten()
    
    for ax_idx, (ax, cov) in enumerate(zip(axes, coverages)):
        pop_data = {}
        all_y_values = []
        
        # Load data for all populations
        for pop in populations:
            file_path = f"{base_path}/header_filt_chr{chr_name}_{pop}_{cov}"
            
            # Check if file exists
            if not os.path.exists(file_path):
                print(f"File not found: {file_path}")
                continue
            
            try:
                # Read the file
                df = pd.read_csv(file_path, sep=r"\s+", header=None,
                                names=["X", "Y"], engine="python")
                df = df.apply(pd.to_numeric, errors="coerce").dropna()
                
                if df.empty:
                    continue
                
                # Process data
                df = detect_chromosome_arms(df)
                df = fix_orientation_and_slope(df)
                
                pop_data[pop] = df
                all_y_values.extend(df["Y"].tolist())
            
            except Exception as e:
                print(f"Error processing {pop}: {e}")
                continue
        
        # Plot data if available
        if pop_data and all_y_values:
            y_min, y_max = min(all_y_values), max(all_y_values)
            
            for pop in populations:
                if pop not in pop_data:
                    continue
                
                df = offset_data(pop_data[pop], pop, cov, y_min, y_max)
                
                # Special styling for 12x coverage
                if cov == "12x":
                    ax.scatter(df["X"], df["Y_offset"],
                             facecolors="none", edgecolors=colors[pop],
                             marker="D", s=marker_size_12x,
                             linewidths=linewidth_12x, zorder=4,
                             alpha=0.9, label=pop if ax_idx == 0 else "")
                else:
                    ax.scatter(df["X"], df["Y_offset"],
                             facecolors=colors[pop], edgecolors="none",
                             marker="o", s=marker_size_small,
                             alpha=0.8, zorder=2,
                             label=pop if ax_idx == 0 else "")
        
        # Subplot formatting
        ax.set_title(f"chr{chr_name} - {cov}", fontsize=12)
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Create legend
    legend_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=colors[pop], markersize=8,
               markeredgecolor='none', label=pop)
        for pop in populations
    ]
    
    # Adjust layout
    plt.subplots_adjust(left=0.08, bottom=0.1, right=0.82, top=0.92,
                       wspace=0.2, hspace=0.3)
    
    # Add legend and axis labels
    fig.legend(handles=legend_handles, loc='upper right',
              bbox_to_anchor=(0.98, 0.55), fontsize=12,
              frameon=True, facecolor='white', edgecolor='black')
    
    fig.text(0.46, 0.02, "Markers", ha="center", fontsize=12)
    fig.text(0.02, 0.5, "Genomes", va="center", rotation="vertical", fontsize=12)
    
    # Save figure
    output_file = f"{output_path}/chr{chr_name}_quadraplot_300dpi.png"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Output for chr{chr_name} saved at {output_file}")

# Plot all chromosomes
def plot_all_chromosomes(chromosomes=None):
    if chromosomes is None:
        chromosomes = [str(i) for i in range(1, 23)] + ['X']
    
    for chr_name in chromosomes:
        try:
            plot_chromosome_quadraplot(chr_name)
        except Exception as e:
            print(f"Error processing chromosome {chr_name}: {e}")

# Main execution
if __name__ == "__main__":
    # Plot chromosome N
    plot_chromosome_quadraplot("N")

