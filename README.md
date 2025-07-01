# ADMIXTURE Plotter

## Overview

ADMIXTURE Plotter is a simple yet effective tool designed for post-processing and visualization of ADMIXTURE results. While we recommend using [pong](https://github.com/ramachandran-lab/pong) for comprehensive ADMIXTURE result visualization, pong does not provide access to the clustered and sorted data files used in its visualizations. This limitation makes it challenging for researchers who need to create customized plots with specific formatting requirements.
ADMIXTURE Plotter addresses this gap by providing a straightforward solution that:
- Aligns ancestry components across different K values based on similarity patterns in Q matrices
- Sorts populations and samples according to their predicted ancestry compositions
- Outputs both adjusted demonstration plots and Q matrices, which enables users to perform customized plotting

## Features

- **Component Alignment**: Automatically aligns ancestry components across different K values using similarity-based algorithms
- **Intelligent Sorting**: Orders populations and samples based on ancestry composition similarity
- **Flexible Input**: Supports multiple Q matrix files and population information formats
- **Customizable Output**: Generates both visualization plots and processed data files for further analysis
- **Population Hierarchy**: Handles two-level population classification (Pop1 and Pop2)

## Usage Guidelines

### 1. Prepare Input Files

Users need to prepare the following files:

- **Q matrices**: ADMIXTURE output files (*.Q) for different K values
- **PLINK bfile**: Specifically the .fam file to determine sample names and order
- **Sample information file** (optional): A TSV format file containing population assignments

### 2. Sample Information File Format

The sample information file should be in TSV format and include the following columns:
- `Sample`: Sample identifier (must match those in the .fam file)
- `Pop1`: Fine-scale population classification
- `Pop2`: Broad-scale population classification

If you only have one level of population classification, set Pop1 and Pop2 to the same values, or provide this information directly in the .fam file.

### 3. Population Ordering

The script will automatically determine optimal population ordering based on ancestry similarity. However, you can also specify custom ordering:

- The script prints the automatically determined Pop1 and Pop2 orders to the console
- These orders are typically not optimal and can be manually adjusted
- Use the printed results as input for `--pop1_order` or `--pop2_order` parameters to achieve desired sample arrangements

### 4. Command Line Usage

```bash
python admixture_plotter.py \
  --fam samples.fam \
  --q_files "prefix.*.Q" \
  --out_prefix output \
  --sample_info sample_info.tsv \
  --pop1_order "Pop1A,Pop1B,Pop1C" \
  --pop2_order "Pop2X,Pop2Y"
  ```

### 5. Output Files

The script generates:

- Visualization plots : PNG and PDF format figures showing aligned ancestry compositions
- Processed Q matrices : TSV files with aligned components and sorted samples
- Console output : Suggested population orders for fine-tuning

## Parameters

- --fam : Path to PLINK .fam file (required)
- --q_files : List of Q matrix files, supports wildcards (required)
- --out_prefix : Output file prefix (required)
- --sample_info : Sample information file in TSV format (optional)
- --pop1_order : Comma-separated Pop1 order (optional)
- --pop2_order : Comma-separated Pop2 order (optional)
- --fig_width : Figure width (default: 20)
- --fig_height_per_k : Figure height per K value (default: 4)

## Example Workflow
1. Run the script with automatic ordering to see suggested population arrangements
2. Adjust the Pop1 and Pop2 orders based on your research needs
3. Re-run the script with custom --pop1_order and --pop2_order parameters
4. Use the output Q matrices for customized plotting and further analysis
