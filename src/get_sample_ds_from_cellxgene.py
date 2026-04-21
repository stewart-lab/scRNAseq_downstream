# %% [markdown]
# # Retrieve a sample of a CellxGene dataset

# %% [markdown]
# This script retrieves and saves a dataset containing a subset of cells from a CellxGene dataset, sampling each annotated cell type. This can serve as the input for `cellxgene_scvi` and other annotation tools to test that they work as expected. It could serve either as a test or a reference dataset. Since CellxGene provides cell type annotations, we have a “gold standard” answer.

# %% [markdown]
# ## Set up

# %% [markdown]
# ### Load libraries

# %%
# Standard python libraries
import numpy as np
import os
import pandas as pd
import warnings

# scVerse
import anndata as ad
from anndata import ImplicitModificationWarning

# CellxGene
import cellxgene_census
import cellxgene_census.experimental

# Local utilities
from utils import load_config, initialize_output_directory, save_package_versions

# Some settings to avoid errors/warnings
# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

# %% [markdown]
# ## Specialized functions for this module

# %% Subsample an AnnData object to get a specified number of cells per cell type
def subsample_by_cell_type(adata, num_cells_per_cell_type):
    """
    Subsample an AnnData object to get a specified number of cells per cell type.
    
    Parameters
    ----------
    adata : AnnData
        The AnnData object to subsample.
    num_cells_per_cell_type : int
        The target number of cells per cell type.
    
    Returns
    -------
    AnnData
        A subsampled AnnData object.
    """
    # Get all unique cell types
    cell_types = adata.obs.loc[adata.obs['cell_type'] != "unknown", 'cell_type'].dropna().unique()
    
    # For each cell type, sample up to num_cells_per_cell_type cells (if available)
    indices_q = []
    
    for ct in cell_types:
        idx = adata.obs[adata.obs['cell_type'] == ct].index.values
        n = min(num_cells_per_cell_type, len(idx))
        if n == 0:
            continue
        # Shuffle indices
        idx = np.random.permutation(idx)
        # If more than num_cells_per_cell_type, select a random set of num_cells_per_cell_type
        if len(idx) >= num_cells_per_cell_type:
            indices_q.extend(list(idx[:num_cells_per_cell_type]))
        else:
            # If less than num_cells_per_cell_type, use all cells for this cell type
            indices_q.extend(list(idx))
    
    # Create the query AnnData subset
    return adata[indices_q, :].copy()

# %% [markdown]
# # Main script
# %% Define the main function
def main():
    # ## Parse the config file and set parameters
    # Load configuration from config.json file
    config_dict = load_config()

    # Initialie output directory
    out_dir = initialize_output_directory(config_dict)

    # Set custom parameters for the method
    # Method name
    METHOD = config_dict["METHOD"]

    # Reference dataset(s)
    census_version = (
        config_dict[METHOD]["reference_datasets"]["census_version"]
    )
    organism = (
        config_dict[METHOD]["reference_datasets"]["organism"]
    )
    ref_dataset_ids = (
        config_dict[METHOD]["reference_datasets"]["ref_dataset_ids"]
    )

    # Number of cells per cell type in a subset of the reference
    ref_cells_per_cell_type = (
        config_dict[METHOD]["ref_cells_per_cell_type"]
    )

    # Name of the output file that would contain the example subset
    output_file = config_dict[METHOD]["output_file"]

    # Random seed
    random_seed = config_dict[METHOD]["random_seed"]

    # ## Load the reference dataset from CellxGene

    # Create a census object
    print("loading census version ", census_version, " ...")
    census = cellxgene_census.open_soma(census_version=census_version)

    # Load the reference dataset
    print("loading reference dataset from CellxGene...")
    adata_census = cellxgene_census.get_anndata(
        census=census,
        measurement_name="RNA",
        organism=organism,
        obs_value_filter=f"dataset_id in {ref_dataset_ids}",
    )
    print(adata_census)

    # ## Sub-sample the data
    # Create a random subset having ```ref_cells_per_cell_type``` representatives of each cell type. 
    print("subsampling the reference dataset...")

    # Set random seed for reproducibility
    np.random.seed(random_seed)

    # Subsample the reference dataset
    adata_subset = subsample_by_cell_type(adata_census, ref_cells_per_cell_type)
    print(adata_subset)

    # ## Save the generated dataset and package versions

    # Save the generated dataset
    print(f"Saving generated dataset to {output_file} ...")
    adata_subset.write_h5ad(out_dir + output_file)

    # Save package versions
    save_package_versions(out_dir)

    # Set permissions for the output directory to be readable and writable by all users
    os.chmod(out_dir, 0o777)

    # All done
    print("Done.")

# %% Run the main function
if __name__ == "__main__":
    main()