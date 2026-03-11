""" 
Utility functions
"""
from datetime import datetime
import json
import matplotlib.pyplot as plt
import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import shutil
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import subprocess

# %% 
# ## Utility functions
# %% Load configuration from config.json file
def load_config():
    """
    Load configuration from config.json file.
    
    Returns
    -------
    dict
        Configuration dictionary loaded from config.json
    """
    work_dir = os.getcwd()
    with open(work_dir + '/config.json') as f:
        config_dict = json.load(f)
        method = config_dict["METHOD"]
        print(
            "loaded parameters from config file: ", 
            config_dict[method]
        )
    return config_dict

# %% Set input data directory based on config settings
# %%
def get_data_dir(config_dict):
    """
    Get the data directory path from configuration.
    
    Parameters
    ----------
    config_dict : dict
        Configuration dictionary containing docker and DATA_DIR settings
    
    Returns
    -------
    str
        Path to the data directory
    """
    docker = config_dict["docker"]
    method = config_dict["METHOD"]
    if docker == "TRUE" or docker == "true" or docker == "T" or docker == "t":
        data_dir = "./data/input_data/"
    else:
        data_dir = config_dict[method]["DATA_DIR"]
    
    return data_dir

# %% Initialize the output directory based on the method name and current timestamp.
def initialize_output_directory(config_dict):
    """
    Initialize the output directory based on the method name and current timestamp.
    
    Parameters
    ----------
    config_dict : dict
        Configuration dictionary
    
    Returns
    -------
    str
        Path to the initialized output directory
    """
    now = datetime.now()
    now = now.strftime("%Y%m%d_%H%M%S")
    out_dir = "./shared_volume/" + config_dict["METHOD"] + "_" + now +"/"
    print("out_dir: ", out_dir)
    os.makedirs(out_dir, mode=0o777, exist_ok=True)
    # copy config file
    shutil.copy('./config.json', out_dir) 

    return out_dir

# %% Save package versions to the output directory
def save_package_versions(out_dir):
    """
    Save package versions from conda and pip to the output directory.
    
    Parameters
    ----------
    out_dir : str
        The output directory where version files will be saved.
    """
    print("Saving package versions...")
    # Check if we're in a conda environment
    in_conda = os.environ.get('CONDA_DEFAULT_ENV') is not None
    if in_conda:
        # Use conda list command
        result = subprocess.run(['conda', 'list', '--explicit'], capture_output=True, text=True)
        packages = result.stdout
        # Write the packages to a file
        with open(out_dir + 'conda-requirements.txt', 'w') as f:
            f.write(packages)
        print("- conda package versions have been written to conda-requirements.txt")

    # Also save pip freeze output
    result = subprocess.run(['pip', 'freeze'], capture_output=True, text=True)
    packages = result.stdout
    # Write the packages to a file
    with open(out_dir + 'pip-requirements.txt', 'w') as f:
        f.write(packages)
    print("- pip package versions have been written to pip-requirements.txt")

# %% Compare two cell metadata columns
def compare_cell_metadata_cols(metadata_col1, metadata_col2, adata, out_dir):
    """
    Compare two cell metadata columns: compute similarity metrics and plot the contingency table.
    
    Parameters
    ----------
    metadata_col1 : str
        The name of the first cell metadata column to compare.
    metadata_col2 : str
        The name of the second cell metadata column to compare.
    adata : AnnData
        Annotated data matrix.
    out_dir : str
        Output directory to save the contingency table plot.
    
    Returns
    -------
    tuple[float, float]
        A tuple containing:
        - ari: Adjusted Rand Index
        - nmi: Normalized Mutual Information
    """

    # Extract the two metadata columns
    col1 = adata.obs[metadata_col1]
    col2 = adata.obs[metadata_col2]

    # Compute and print similarity metrics
    ari = adjusted_rand_score(col1, col2)
    nmi = normalized_mutual_info_score(col1, col2)
    print()
    print(f"Comparing '{metadata_col1}' and '{metadata_col2}':")
    print(f"Adjusted Rand Index (ARI): {ari:.4f}")
    print(f"Normalized Mutual Information (NMI): {nmi:.4f}")

    # Create a contingency table
    contingency_table = pd.crosstab(col1, col2)

    # Decide how to annotate the heatmap based on table size to avoid clutter
    n_rows, n_cols = contingency_table.shape
    max_dim = max(n_rows, n_cols)
    n_cells = n_rows * n_cols

    # Default: annotate with a readable font size for small tables
    annot = True
    heatmap_kwargs = {}

    # For very large tables, disable annotations entirely
    if max_dim > 50 or n_cells > 1000:
        annot = False
    # For moderately large tables, keep annotations but shrink font size
    elif max_dim > 20 or n_cells > 400:
        heatmap_kwargs["annot_kws"] = {"size": 6}
    else:
        heatmap_kwargs["annot_kws"] = {"size": 8}

    # Plot the contingency table
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        contingency_table,
        annot=annot,
        fmt='d' if annot else '',
        cmap='viridis',
        **heatmap_kwargs,
    )
    plt.title(f'Contingency Table: {metadata_col1} vs {metadata_col2}\nARI: {ari:.4f}, NMI: {nmi:.4f}')
    plt.xlabel(metadata_col2)
    plt.ylabel(metadata_col1)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(out_dir, f'{metadata_col1}_vs_{metadata_col2}_contingency.png'), bbox_inches='tight')
    plt.close()

    return ari, nmi

# %% Plot UMAPs with consistent figure dimensions and spacing

# Default layout parameters for UMAP plots (in inches and dpi)
UMAP_LAYOUT_DEFAULTS = {
    "fixed_dpi": 300,
    "plot_area_in": 5.2,
    "left_margin_in": 0.55,
    "bottom_margin_in": 0.55,
    "top_margin_in": 0.35,
    "right_margin_continuous_in": 1.5,
    "right_margin_categorical_in": 4.0,
}

# The UMAP plotting function
# Example call: override one default layout value, keep all others at default
# plot_umaps(adata, out_dir, plot_specs, layout={"right_margin_categorical_in": 3.0})
def plot_umaps(adata, out_dir, plot_specs, layout=None):
    """
    Plot UMAPs with consistent figure dimensions and spacing.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data to plot.
    out_dir : str
        Directory to save the output plots.
    plot_specs : list[tuple[str, str, str]]
        List of (color_col, filename, title) tuples specifying what to plot.
    layout : dict | None, default None
        Optional overrides for UMAP_LAYOUT_DEFAULTS.
        Example: {"right_margin_categorical_in": 3.0}
    """
    cfg = UMAP_LAYOUT_DEFAULTS | (layout or {})

    for color_col, filename, title in plot_specs:
        series = adata.obs[color_col]
        is_continuous = pd.api.types.is_numeric_dtype(series)

        right_margin_in = (
            cfg["right_margin_continuous_in"]
            if is_continuous
            else cfg["right_margin_categorical_in"]
        )

        fig_w = cfg["left_margin_in"] + cfg["plot_area_in"] + right_margin_in
        fig_h = cfg["bottom_margin_in"] + cfg["plot_area_in"] + cfg["top_margin_in"]
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=cfg["fixed_dpi"])

        ax_left = cfg["left_margin_in"] / fig_w
        ax_bottom = cfg["bottom_margin_in"] / fig_h
        ax_width = cfg["plot_area_in"] / fig_w
        ax_height = cfg["plot_area_in"] / fig_h

        with plt.rc_context():
            sc.pl.umap(
                adata,
                color=color_col,
                title=title,
                ax=ax,
                show=False,
                legend_loc="right margin",
                colorbar_loc="right",
            )

            ax.set_position([ax_left, ax_bottom, ax_width, ax_height])
            ax.set_box_aspect(1)

            fig.savefig(
                os.path.join(out_dir, filename),
                dpi=cfg["fixed_dpi"],
                bbox_inches="tight",
                pad_inches=0.02,
            )
            plt.close(fig)
