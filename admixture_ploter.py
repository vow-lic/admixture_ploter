import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import glob
import argparse
import re

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='ADMIXTURE results visualization and processing tool')
    
    # Required arguments
    parser.add_argument('--fam', required=True, help='Path to PLINK FAM file')
    parser.add_argument('--q_files', required=True, nargs='+', help='Q matrix files list, supports multiple files or wildcard patterns (e.g., "prefix.5.Q prefix.6.Q" or "prefix.*.Q")')
    parser.add_argument('--out_prefix', required=True, help='Output file prefix')
    
    # Optional arguments
    parser.add_argument('--sample_info', help='Sample information file path (TSV format with Sample, Pop1, Pop2 columns)')
    parser.add_argument('--pop1_order', help='Pop1 order, comma-separated, fine-scale population classifications')
    parser.add_argument('--pop2_order', help='Pop2 order, comma-separated, broad-scale population classifications')
    parser.add_argument('--fig_width', type=float, default=20, help='Figure width in inches')
    parser.add_argument('--fig_height_per_k', type=float, default=4, help='Figure height per K value in inches')
    
    return parser.parse_args()

def read_files(args):
    """Read input files"""
    # Read Q files
    q_files = {}
    k_list = []
    
    # Handle wildcards in q_files argument
    q_file_paths = []
    for q_file_arg in args.q_files:
        # If argument contains wildcard, expand it
        if '*' in q_file_arg:
            expanded_files = glob.glob(q_file_arg)
            q_file_paths.extend(expanded_files)
        else:
            # Otherwise treat as a direct file path
            q_file_paths.append(q_file_arg)
    
    print(f"Found {len(q_file_paths)} Q files")
    
    # Read fam file to get sample IDs
    fam_data = pd.read_csv(args.fam, sep=r'\s+', header=None)
    sample_ids = fam_data.iloc[:, 1].tolist()
    print(f"Read fam file: {len(sample_ids)} sample IDs")
    
    # Extract K values from Q file names and read Q files
    for q_file_path in q_file_paths:
        if not os.path.exists(q_file_path):
            print(f"Warning: Q file not found: {q_file_path}")
            continue
        
        # Extract K value from filename using regex
        k_match = re.search(r'\.(\d+)\.Q$', q_file_path)
        if k_match:
            k = int(k_match.group(1))
            k_list.append(k)
            
            # Read Q file
            q_data = pd.read_csv(q_file_path, sep=r'\s+', header=None)
            q_data.columns = [f"K{i+1}" for i in range(q_data.shape[1])]
            q_data['sample_id'] = sample_ids
            q_files[k] = q_data
            print(f"Read K={k} Q file: {q_data.shape[1]-1} ancestry components")
        else:
            print(f"Warning: Could not extract K value from filename: {q_file_path}")
    
    # Sort k_list
    k_list.sort()
    
    # Read sample information if provided
    if args.sample_info:
        pop_info = pd.read_csv(args.sample_info, sep='\t')
        print(f"Read sample information: {len(pop_info)} samples")
        
        # Merge sample information with Q files
        for k in k_list:
            q_files[k] = pd.merge(q_files[k], pop_info, left_on='sample_id', right_on='Sample', how='left')
            # Fill missing Pop1 and Pop2 information with " "
            q_files[k]['Pop1'] = q_files[k]['Pop1'].fillna(' ')
            q_files[k]['Pop2'] = q_files[k]['Pop2'].fillna(' ')
    else:
        # If no sample info provided, use fam file's first column as Pop1 and Pop2
        print("No sample information provided, using fam file's first column as Pop1 and Pop2")
        for k in k_list:
            q_files[k]['Pop1'] = fam_data.iloc[:, 0].values
            q_files[k]['Pop2'] = fam_data.iloc[:, 0].values
    
    # Parse Pop1 and Pop2 order
    pop1_order = args.pop1_order.split(',') if args.pop1_order else []
    pop2_order = args.pop2_order.split(',') if args.pop2_order else []
    
    return k_list, q_files, pop1_order, pop2_order

def order_samples(k_list, q_files, pop1_order, pop2_order):
    """Order samples based on ancestry components"""
    # Use the minimum K value result to determine ordering
    min_k = min(k_list)
    min_q_data = q_files[min_k]
    
    # Get all unique Pop1 and Pop2 values from the data
    all_pop1 = set(q_files[max(k_list)]['Pop1'].unique())
    all_pop2 = set(min_q_data['Pop2'].unique())
    
    # Remove space placeholder if present
    if ' ' in all_pop1:
        all_pop1.remove(' ')
    if ' ' in all_pop2:
        all_pop2.remove(' ')
    
    # Check if provided Pop2 order is valid
    if pop2_order:
        # Check if all provided Pop2 values exist in the data
        invalid_pop2 = set(pop2_order) - all_pop2
        if invalid_pop2:
            print(f"Warning: The following Pop2 values in pop2_order do not exist in the data: {', '.join(invalid_pop2)}")
        
        # Check if all Pop2 values in the data are included in the provided order
        missing_pop2 = all_pop2 - set(pop2_order)
        if missing_pop2:
            print(f"Error: The following Pop2 values from the data are missing in pop2_order: {', '.join(missing_pop2)}")
            print("Please include all Pop2 values in pop2_order or leave it empty for automatic ordering.")
            exit(1)
    
    # Check if provided Pop1 order is valid
    if pop1_order:
        # Check if all provided Pop1 values exist in the data
        invalid_pop1 = set(pop1_order) - all_pop1
        if invalid_pop1:
            print(f"Warning: The following Pop1 values in pop1_order do not exist in the data: {', '.join(invalid_pop1)}")
        
        # Check if all Pop1 values in the data are included in the provided order
        missing_pop1 = all_pop1 - set(pop1_order)
        if missing_pop1:
            print(f"Error: The following Pop1 values from the data are missing in pop1_order: {', '.join(missing_pop1)}")
            print("Please include all Pop1 values in pop1_order or leave it empty for automatic ordering.")
            exit(1)

    # Calculate the main component for each Pop2
    pop2_component = {}
    for pop2_name, group in min_q_data.groupby('Pop2'):
        if pop2_name == ' ':
            continue
        # Select K columns (only numerical columns)
        k_columns = [col for col in group.columns if col.startswith('K')]
        # Calculate average ancestry proportions for this Pop2
        avg_proportions = group[k_columns].mean().values
        pop2_component[pop2_name] = avg_proportions
    
    # Find the K value with highest proportion across all samples
    k_columns = [col for col in min_q_data.columns if col.startswith('K')]
    all_proportions = min_q_data[k_columns].mean().values
    kmax = np.argmax(all_proportions)
    
    # Order Pop2
    if pop2_order:  # If Pop2 order is provided
        pop2_order_final = pop2_order
        print(f"Using preset Pop2 order: {pop2_order_final}")
    else:
        pop2_order_final = []
        remaining_pop2s = list(pop2_component.keys())
        
        # 1. Find Pop2 with highest KMAX proportion as the first
        if remaining_pop2s:
            first_pop2 = max(remaining_pop2s, key=lambda x: pop2_component[x][kmax])
            pop2_order_final.append(first_pop2)
            remaining_pop2s.remove(first_pop2)
        
        # 2. Iteratively find the most similar Pop2 to the last one
        while remaining_pop2s:
            last_pop2 = pop2_order_final[-1]
            # Calculate similarity to the last Pop2 (Euclidean distance)
            similarities = [
                (pop2, np.linalg.norm(pop2_component[pop2] - pop2_component[last_pop2]))
                for pop2 in remaining_pop2s
            ]
            # Find the most similar Pop2
            next_pop2 = min(similarities, key=lambda x: x[1])[0]
            pop2_order_final.append(next_pop2)
            remaining_pop2s.remove(next_pop2)
        
        # Add unknown category
        if ' ' in min_q_data['Pop2'].unique():
            pop2_order_final.append(' ')
        
        print(f"Pop2 order: {','.join(pop2_order_final)}")
    
    # Order Pop1 and samples based on maximum K value
    max_k = max(k_list)
    max_q_data = q_files[max_k]
    
    # Calculate the main component for each Pop1
    pop1_component = {}
    for pop1_name, group in max_q_data.groupby('Pop1'):
        if pop1_name == ' ':
            continue
        # Select K columns (only numerical columns)
        k_columns = [col for col in group.columns if col.startswith('K')]
        # Calculate average ancestry proportions for this Pop1
        avg_proportions = group[k_columns].mean().values
        pop1_component[pop1_name] = avg_proportions
    
    # Order Pop1 by Pop2 and similarity
    if pop1_order:  # If Pop1 order is provided
        pop1_order_final = pop1_order
        print(f"Using preset Pop1 order: {pop1_order_final}")
    else:
        pop1_order_final = []
        for pop2_name in pop2_order_final:
            if pop2_name == ' ':
                continue
                
            # Get all Pop1s under current Pop2
            pop2_data = max_q_data[max_q_data['Pop2'] == pop2_name]
            pop2_pop1s = list(pop2_data['Pop1'].unique())
            if ' ' in pop2_pop1s:
                pop2_pop1s.remove(' ')
            
            if not pop2_pop1s:
                continue
                
            # Calculate average ancestry proportions for each Pop1 in current Pop2
            pop2_pop1_proportions = {}
            for pop1 in pop2_pop1s:
                pop1_data = pop2_data[pop2_data['Pop1'] == pop1]
                k_columns = [col for col in pop1_data.columns if col.startswith('K')]
                pop2_pop1_proportions[pop1] = pop1_data[k_columns].mean().values
            
            # Order Pop1s within current Pop2
            pop2_pop1_order = []
            remaining_pop1s = pop2_pop1s.copy()
            
            # 1. Find Pop1 with highest KMAX proportion as the first
            if remaining_pop1s:
                first_pop1 = max(remaining_pop1s, key=lambda x: pop2_pop1_proportions[x][kmax])
                pop2_pop1_order.append(first_pop1)
                remaining_pop1s.remove(first_pop1)
            
            # 2. Iteratively find the most similar Pop1 to the last one
            while remaining_pop1s:
                last_pop1 = pop2_pop1_order[-1]
                # Calculate similarity to the last Pop1 (Euclidean distance)
                similarities = [
                    (pop1, np.linalg.norm(pop2_pop1_proportions[pop1] - pop2_pop1_proportions[last_pop1]))
                    for pop1 in remaining_pop1s
                ]
                # Find the most similar Pop1
                next_pop1 = min(similarities, key=lambda x: x[1])[0]
                pop2_pop1_order.append(next_pop1)
                remaining_pop1s.remove(next_pop1)
            
            # Add current Pop2's Pop1 order to the total order
            pop1_order_final.extend(pop2_pop1_order)
        
        # Add unknown category
        if ' ' in max_q_data['Pop1'].unique():
            pop1_order_final.append(' ')
        
        print(f"Pop1 order: {','.join(pop1_order_final)}")
    
    # Order samples by Pop2, Pop1, and individual ancestry proportions
    sample_order = []
    for pop2_name in pop2_order_final:
        pop2_data = max_q_data[max_q_data['Pop2'] == pop2_name]
        for pop1_name in pop1_order_final:
            pop1_data = pop2_data[pop2_data['Pop1'] == pop1_name]
            if pop1_data.empty:
                continue
            
            # Get ancestry proportions for all samples in current Pop1
            k_columns = [col for col in pop1_data.columns if col.startswith('K')]
            sample_proportions = {}
            for _, row in pop1_data.iterrows():
                sample_proportions[row['sample_id']] = row[k_columns].values
            
            # Order samples within current Pop1
            pop1_sample_order = []
            remaining_samples = list(sample_proportions.keys())
            
            # 1. Find sample with highest KMAX proportion as the first
            if remaining_samples:
                first_sample = max(remaining_samples, 
                                 key=lambda x: sample_proportions[x][kmax])
                pop1_sample_order.append(first_sample)
                remaining_samples.remove(first_sample)
            
            # 2. Iteratively find the most similar sample to the last one
            while remaining_samples:
                last_sample = pop1_sample_order[-1]
                # Calculate similarity to the last sample (Euclidean distance)
                similarities = [
                    (sample, np.linalg.norm(sample_proportions[sample] - 
                                          sample_proportions[last_sample]))
                    for sample in remaining_samples
                ]
                # Find the most similar sample
                next_sample = min(similarities, key=lambda x: x[1])[0]
                pop1_sample_order.append(next_sample)
                remaining_samples.remove(next_sample)
            
            # Add current Pop1's sample order to the total order
            sample_order.extend(pop1_sample_order)
    
    # Add unknown category samples
    unknown_samples = max_q_data[~max_q_data['sample_id'].isin(sample_order)]['sample_id'].tolist()
    sample_order.extend(unknown_samples)
    
    print(f"Sample ordering completed: {len(sample_order)} samples")
    
    return pop2_order_final, pop1_order_final, sample_order, kmax

def map_components(k_list, q_files):
    """Map components across different K values"""
    def calculate_similarity(dist1, dist2):
        """Calculate similarity between two distributions"""
        return np.linalg.norm(dist1 - dist2)
    
    # Process color mapping
    ck_mapping = {}  # Store K to CK mapping for each K value file
    
    # 1. Process the maximum K value file
    max_k = max(k_list)
    max_k_data = q_files[max_k]
    max_k_columns = [col for col in max_k_data.columns if col.startswith('K')]
    
    # Create CK mapping for maximum K value
    ck_mapping[max_k] = {f"K{i+1}": f"CK{i+1}" for i in range(max_k)}
    # Calculate component distributions for maximum K value
    ck_distributions = {}
    for k_col in max_k_columns:
        ck = ck_mapping[max_k][k_col]
        ck_distributions[ck] = max_k_data[k_col].values
    
    # 2. Process other K values, from large to small
    for k in sorted(k_list, reverse=True):
        if k == max_k:
            continue
            
        q_data = q_files[k]
        k_columns = [col for col in q_data.columns if col.startswith('K')]
        ck_mapping[k] = {}
        
        # Calculate component distributions for current K value
        k_distributions = {col: q_data[col].values for col in k_columns}
        
        # Find the best matching CK for each K
        used_cks = set()
        for k_col in k_columns:
            k_dist = k_distributions[k_col]
            
            # Calculate similarity to all available CKs
            similarities = []
            for ck, ck_dist in ck_distributions.items():
                if ck not in used_cks:
                    sim = calculate_similarity(k_dist, ck_dist)
                    similarities.append((ck, sim))
            
            # Find the best matching CK
            if similarities:
                best_ck = min(similarities, key=lambda x: x[1])[0]
                ck_mapping[k][k_col] = best_ck
                used_cks.add(best_ck)
                
                # Update CK distribution (take average)
                ck_distributions[best_ck] = (ck_distributions[best_ck] + k_dist) / 2
            else:
                # If no available CK, create a new one
                new_ck = f"CK{len(ck_distributions) + 1}"
                ck_mapping[k][k_col] = new_ck
                ck_distributions[new_ck] = k_dist
    
    # Print mapping results for each K value
    for k in sorted(k_list):
        print(f"\nMapping results for K={k}:")
        mapping = ck_mapping[k]
        sorted_items = sorted(mapping.items(), key=lambda x: int(x[0][1:]))
        for k_col, ck in sorted_items:
            print(f"{k_col} -> {ck}")
    
    return ck_mapping

def plot_admixture(k_list, q_files, pop2_order, pop1_order, sample_order, ck_mapping, output_prefix, fig_width, fig_height_per_k):
    """Plot ADMIXTURE results"""
    # Create color scheme
    max_ck = max([int(ck[2:]) for ck in set().union(*[set(mapping.values()) for mapping in ck_mapping.values()])])
    color_palette = sns.color_palette("husl", max_ck)
    ck_colors = {f"CK{i+1}": color_palette[i] for i in range(max_ck)}
    
    # Create figure
    fig = plt.figure(figsize=(fig_width, fig_height_per_k * len(k_list)))
    gs = gridspec.GridSpec(len(k_list), 1, height_ratios=[1] * len(k_list))
    
    # Record starting positions for each Pop1 and Pop2, for annotation
    pop1_positions = {}
    pop2_positions = {}
    
    # Create a subplot for each K value
    for i, k in enumerate(sorted(k_list)):
        if k not in q_files:
            continue
            
        ax = plt.subplot(gs[i])
        
        # Get data for current K value
        q_data = q_files[k]
        
        # Reorder data according to sorted sample order
        q_data_sorted = q_data.set_index('sample_id').loc[sample_order].reset_index()
        
        # Draw stacked bar chart
        bottom = np.zeros(len(sample_order))
        k_columns = [col for col in q_data.columns if col.startswith('K')]
        
        for k_col in k_columns:
            # Use color of mapped CK
            ck = ck_mapping[k][k_col]
            color = ck_colors[ck]
            ax.bar(range(len(sample_order)), q_data_sorted[k_col], 
                   bottom=bottom, width=1, color=color)
            bottom += q_data_sorted[k_col].values
        
        # Set y-axis range and label
        ax.set_ylim(0, 1)
        ax.set_ylabel(f"K={k}", fontsize=14)
        
        # Remove x-axis ticks and labels
        ax.set_xticks([])
        ax.set_xticklabels([])
        
        # Add Pop2 dividing lines
        current_pop2 = None
        for idx, sample in enumerate(sample_order):
            sample_info = q_data[q_data['sample_id'] == sample]
            pop2 = sample_info['Pop2'].iloc[0] if not sample_info.empty else ' '
            
            if pop2 != current_pop2:
                ax.axvline(x=idx-0.5, color='black', linestyle='-', 
                          linewidth=0.8, alpha=0.7)
                current_pop2 = pop2
        
        # Add Pop1 dividing lines (thinner and lighter than Pop2 lines)
        current_pop1 = None
        for idx, sample in enumerate(sample_order):
            sample_info = q_data[q_data['sample_id'] == sample]
            pop1 = sample_info['Pop1'].iloc[0] if not sample_info.empty else ' '
            
            if pop1 != current_pop1:
                ax.axvline(x=idx-0.5, color='gray', linestyle='-', 
                          linewidth=0.3, alpha=0.5)
                current_pop1 = pop1
                
        # Only record Pop1 and Pop2 positions in the first iteration
        if i == 0:
            current_pop2 = None
            current_pop1 = None
            for idx, sample in enumerate(sample_order):
                sample_info = q_data[q_data['sample_id'] == sample]
                pop1 = sample_info['Pop1'].iloc[0] if not sample_info.empty else ' '
                pop2 = sample_info['Pop2'].iloc[0] if not sample_info.empty else ' '
                
                # Record Pop1 starting position
                if pop1 != current_pop1:
                    pop1_positions[pop1] = idx
                    current_pop1 = pop1
                
                # Record Pop2 starting position
                if pop2 != current_pop2:
                    pop2_positions[pop2] = idx
                    current_pop2 = pop2
        
        # Save data for this subplot
        output_data = q_data_sorted[['sample_id', 'Pop1', 'Pop2'] + k_columns].copy()
        
        # Rename K columns to CK columns
        for k_col in k_columns:
            ck = ck_mapping[k][k_col]
            output_data.rename(columns={k_col: ck}, inplace=True)
        
        # Save to file
        output_file = f"{output_prefix}.{i+1}.{k}.plotQ.tsv"
        output_data.to_csv(output_file, sep='\t', index=False)
        print(f"Data for subplot {i+1} (K={k}) saved to {output_file}")
    
    # Add Pop1 and Pop2 labels below the last subplot
    ax = plt.subplot(gs[-1])
    
    # Calculate maximum Pop1 label length for dynamic Pop2 positioning
    max_pop1_length = 0
    displayed_pop1s = []
    
    # Add Pop1 labels (only for Pop1s with more than 5 samples)
    for i, pop1 in enumerate(pop1_order):
        if pop1 in pop1_positions:
            next_pop1 = pop1_order[i+1] if i < len(pop1_order)-1 else None
            end_pos = pop1_positions.get(next_pop1, len(sample_order))
            pop1_samples_count = end_pos - pop1_positions[pop1]
            
            if pop1_samples_count >= 5:  # Only label Pop1s with more than 5 samples
                mid_pos = pop1_positions[pop1] + (end_pos - pop1_positions[pop1])/2
                ax.annotate(pop1, xy=(mid_pos, -0.05), xycoords=('data', 'axes fraction'), 
                           ha='center', va='top', fontsize=8, rotation=-90)
                displayed_pop1s.append(pop1)
                max_pop1_length = max(max_pop1_length, len(pop1))
    
    # Calculate dynamic Pop2 label position based on Pop1 label length
    # Base position is -0.15, adjust downward based on Pop1 string length
    # Each character adds approximately 0.016 to the offset (empirically determined)
    pop2_y_offset = -0.15 - (max_pop1_length * 0.016)
    
    # Add Pop2 labels
    for i, pop2 in enumerate(pop2_order):
        if pop2 in pop2_positions:
            next_pop2 = pop2_order[i+1] if i < len(pop2_order)-1 else None
            end_pos = pop2_positions.get(next_pop2, len(sample_order))
            mid_pos = pop2_positions[pop2] + (end_pos - pop2_positions[pop2])/2
            ax.annotate(pop2, xy=(mid_pos, pop2_y_offset), xycoords=('data', 'axes fraction'), 
                       ha='center', va='top', fontsize=12, fontweight='bold')
    
    # Add legend
    handles = []
    labels = []
    for i in range(max_ck):
        ck = f"CK{i+1}"
        if ck in ck_colors:
            handles.append(plt.Rectangle((0,0), 1, 1, color=ck_colors[ck]))
            labels.append(f"Ancestry {i+1}")
    
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.02), 
              ncol=min(10, max_ck), fontsize=10)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    plt.subplots_adjust(hspace=0.1)
    
    # Save figure
    plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.pdf", bbox_inches='tight')
    print(f"Figures saved as {output_prefix}.png and {output_prefix}.pdf")

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Read input files
    k_list, q_files, pop1_order, pop2_order = read_files(args)
    
    # If no Q files were successfully read, exit program
    if not q_files:
        print("Error: No Q files were successfully read, program exiting")
        return
    
    # Order samples
    pop2_order_final, pop1_order_final, sample_order, kmax = order_samples(k_list, q_files, pop1_order, pop2_order)
    
    # Map components across different K values
    ck_mapping = map_components(k_list, q_files)
    
    # Plot ADMIXTURE results
    plot_admixture(k_list, q_files, pop2_order_final, pop1_order_final, sample_order, ck_mapping, 
                  args.out_prefix, args.fig_width, args.fig_height_per_k)

if __name__ == "__main__":
    main()