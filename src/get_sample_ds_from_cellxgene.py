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
from datetime import datetime
import json
import numpy as np
import os
import pandas as pd
import shutil
import subprocess
import warnings

# scVerse
import anndata as ad
from anndata import ImplicitModificationWarning

# CellxGene
import cellxgene_census
import cellxgene_census.experimental

# Some settings to avoid errors/warnings
# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

# %% [markdown]
# ### Parse the config file.

# %%
# Read in the config file
work_dir = os.getcwd()
with open(work_dir+'/config.json') as f:
    config_dict = json.load(f)
    print(
        "loaded config file: ", 
        config_dict["get_sample_ds_from_cellxgene"]
    )

# docker and data directory
# (the data directory is not really used in this script, but we parse it for 
# consistency with other pipeline scripts)
docker = config_dict["docker"]
if docker == "TRUE" or docker == "true" or docker == "T" or docker == "t":
    DATA_DIR = "./data/input_data/"
else:
    DATA_DIR = config_dict["get_sample_ds_from_cellxgene"]["DATA_DIR"]

# Reference dataset(s)
census_version = (
    config_dict["get_sample_ds_from_cellxgene"]["reference_datasets"]["census_version"]
)
organism = (
    config_dict["get_sample_ds_from_cellxgene"]["reference_datasets"]["organism"]
)
ref_dataset_ids = (
    config_dict["get_sample_ds_from_cellxgene"]["reference_datasets"]["ref_dataset_ids"]
)

# Number of cells per cell type in a subset of the reference
ref_cells_per_cell_type = (
    config_dict["get_sample_ds_from_cellxgene"]["ref_cells_per_cell_type"]
)

# Name of the output file that would contain the example subset
output_file = config_dict["get_sample_ds_from_cellxgene"]["output_file"]

# Random seed
random_seed = config_dict["get_sample_ds_from_cellxgene"]["random_seed"]

# %% [markdown]
# ### Initialize the output directory

# %%
print("Initializing the output directory...")
now = datetime.now()
now = now.strftime("%Y%m%d_%H%M%S")
out_dir = "./shared_volume/" + config_dict["METHOD"] + "_" + now +"/"
print("out_dir: ", out_dir)
os.makedirs(out_dir, mode=0o777, exist_ok=True)
# copy config file
shutil.copy('./config.json', out_dir) 

# %% [markdown]
# ### Load the reference dataset from CellxGene

# %% [markdown]
# Create a census object

# %%
print("loading census version ", census_version, " ...")
census = cellxgene_census.open_soma(census_version=census_version)

# %% [markdown]
# Load the reference dataset

# %%
print("loading reference dataset from CellxGene...")
adata_census = cellxgene_census.get_anndata(
    census=census,
    measurement_name="RNA",
    organism=organism,
    obs_value_filter=f"dataset_id in {ref_dataset_ids}",
)
print(adata_census)

# %% [markdown]
# ## Sub-sample the data
# Create a random subset having ```ref_cells_per_cell_type``` representatives of each cell type. 

# %%
print("subsampling the reference dataset...")

def subsample_by_cell_type(adata, num_cells_per_cell_type):
    """
    Subsample an AnnData object to have a specified number of cells per cell type.
    
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

# Set random seed for reproducibility
np.random.seed(random_seed)

# Subsample the reference dataset
adata_query = subsample_by_cell_type(adata_census, ref_cells_per_cell_type)
print(adata_query)

# %% [markdown]
# ## Save the generated dataset and package versions

# %%
# Save the generated dataset
print(f"Saving generated dataset to {output_file} ...")
adata_query.write_h5ad(out_dir + output_file)

# save package versions
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

os.chmod(out_dir, 0o777)

print("Done.")