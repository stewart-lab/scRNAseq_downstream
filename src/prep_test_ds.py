"""
Add noise and modify cell ids in a test dataset
"""

# %% Import libraries
# Standard python libraries
import numpy as np
import os
import pandas as pd
import warnings

# scVerse
import anndata as ad
from anndata import ImplicitModificationWarning
import scanpy as sc
import scipy.sparse as sp

# Local utilities
from utils import load_config, get_data_dir, initialize_output_directory, save_package_versions

# Some settings to avoid errors/warnings
# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

# %% [markdown]
# ## Specialized functions for this module
def add_noise_to_counts(adata, noise_level=0.067):
    """
    Add noise to the count data in an AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        The AnnData object to modify.
    noise_level : float, optional
        The level of noise to add (default is 0.067).
    
    Returns
    -------
    AnnData
        An AnnData object with noise added to the count data.
    """
    # Get the count matrix as a numpy array
    if sp.issparse(adata.X):
        counts = adata.X.toarray()
    else:
        counts = adata.X
    
    # Add noise to the counts
    noise = np.random.normal(
        loc=0.0, 
        scale=noise_level * np.mean(counts, axis=0, keepdims=True), 
        size=counts.shape
    )
    noisy_counts = counts + np.round(noise)
    
    # Ensure that counts remain non-negative and convert back to sparse format if needed
    noisy_counts = np.clip(noisy_counts, a_min=0, a_max=None)
    
    if isinstance(adata.X, np.ndarray):
        adata.X = noisy_counts
    else:
        adata.X = sp.csr_matrix(noisy_counts)
    
    return adata

# %% Modify cell ids by adding a prefix
def modify_cell_ids(adata, prefix="test_ds"):
    """
    Modify cell IDs in an AnnData object by adding a prefix.
    
    Parameters
    ----------
    adata : AnnData
        The AnnData object to modify.
    prefix : str, optional
        The prefix to add to cell IDs (default is "test_ds").
    
    Returns
    -------
    AnnData
        An AnnData object with modified cell IDs.
    """
    adata.obs_names = [f"{prefix}_{cell_id}" for cell_id in adata.obs_names]
    return adata

# %% [markdown]
# # Main script
# %% Define the main function
def main():
    # ## Parse the config file and set parameters
    # Load configuration from config.json file
    config_dict = load_config()

    # Set custom parameters for the method
    # Method name
    METHOD = config_dict["METHOD"]

    # Input file
    DATA_DIR = get_data_dir(config_dict)
    input_file = DATA_DIR + config_dict[METHOD]["input_file"]

    # Output file
    out_dir = initialize_output_directory(config_dict)
    output_file = config_dict[METHOD]["output_file"]

    # Noise level 
    noise_level = config_dict[METHOD]["noise_level"]

    # Random seed
    random_seed = config_dict[METHOD]["random_seed"]

    # ### Load the input dataset
    print(f"Loading input dataset from file: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(adata)

    # Deep copy the adata object to avoid modifying the original data
    adata2 = adata.copy()

    # ## Add noise to the count data
    print("adding noise to the count data...")
    # Set random seed for reproducibility
    np.random.seed(random_seed)
    # Add noise to the count data
    adata2 = add_noise_to_counts(adata2, noise_level=noise_level)
    print(adata2)

    # ## Modify cell ids by adding a prefix
    print("modifying cell ids by adding a prefix...")
    adata2 = modify_cell_ids(adata2, prefix="test_ds")
    print(adata2)

    # ## Visualize the modified dataset compared to the original using UMAP
    print("visualizing the modified dataset compared to the original using UMAP...")

    # Combine the original and modified datasets for visualization
    adata_combined = ad.concat([adata, adata2], axis=0)
    adata_combined.obs["dataset"] = ["original"] * adata.shape[0] + ["modified"] * adata2.shape[0]

    # Normalize and log-transform the combined dataset for visualization
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)

    # Select highly variable genes for visualization
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, batch_key="dataset")

    # Perform PCA and UMAP for visualization
    sc.pp.pca(adata_combined)
    sc.pp.neighbors(adata_combined, use_rep="X_pca")
    sc.tl.umap(adata_combined)
    (sc.pl.umap(
        adata_combined, color="dataset", 
        title="UMAP of Original vs Modified Dataset", show=False).
        figure.savefig(out_dir + "umap_original_vs_modified.png")
    )

    # ## Save the generated dataset and package versions

    # Save the generated dataset
    print(f"Saving generated dataset to {output_file} ...")
    adata2.write_h5ad(out_dir + output_file)

    # Save package versions
    save_package_versions(out_dir)

    # Set permissions for the output directory to be readable and writable by all users
    os.chmod(out_dir, 0o777)

    # All done
    print("Done.")

# %% Run the main function
if __name__ == "__main__":
    main()