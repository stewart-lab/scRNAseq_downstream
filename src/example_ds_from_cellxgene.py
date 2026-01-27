# %% [markdown]
# # Make a small example dataset from CellxGene

# %% [markdown]
# This can serve as the input for `cellxgene_scvi` and other annotation tools to test that they work as expected. This dataset will contain a subset of cells from a reference dataset in CellxGene. Since CellxGene provides cell type annotations, we have a “gold standard” answer.

# %% [markdown]
# ## Set up

# %% [markdown]
# ### Load libraries

# %%
print("loading libraries...")

# Standard python libraries
from datetime import datetime
import json
import numpy as np
import os
import pandas as pd
import shutil
import subprocess

# scVerse
import anndata

# CellxGene
import cellxgene_census
import cellxgene_census.experimental

# Ontology
from oaklib import get_adapter


# %% [markdown]
# Make sure the plots show up in the notebook

# %%
# %matplotlib inline

# %% [markdown]
# ### Parse the config file.

# %%
# Read in the config file
work_dir = os.getcwd()
with open(work_dir+'/config.json') as f:
    config_dict = json.load(f)
    print(
        "loaded config file: ", 
        config_dict["example_ds_from_cellxgene"]
    )

# docker and data directory
# (the data directory is not really used in this script, but we parse it for 
# consistency with other pipeline scripts)
docker = config_dict["docker"]
if docker == "TRUE" or docker == "true" or docker == "T" or docker == "t":
    DATA_DIR = "./data/input_data/"
else:
    DATA_DIR = config_dict["example_ds_from_cellxgene"]["DATA_DIR"]

# Reference dataset(s)
census_version = (
    config_dict["example_ds_from_cellxgene"]["reference_datasets"]["census_version"]
)
organism = (
    config_dict["example_ds_from_cellxgene"]["reference_datasets"]["organism"]
)
ref_dataset_ids = (
    config_dict["example_ds_from_cellxgene"]["reference_datasets"]["ref_dataset_ids"]
)

# High-level cell types
high_level_cell_types = (
    config_dict["example_ds_from_cellxgene"]["high_level_cell_types"]
)

# Number of cells per cell type in a subset of the reference
ref_cells_per_cell_type = (
    config_dict["example_ds_from_cellxgene"]["ref_cells_per_cell_type"]
)

# Name of the output file that would contain the example subset
output_file = config_dict["example_ds_from_cellxgene"]["output_file"]

# Random seed
random_seed = config_dict["example_ds_from_cellxgene"]["random_seed"]

# %% [markdown]
# ### Initialize the output directory

# %%
now = datetime.now()
now = now.strftime("%Y%m%d_%H%M%S")
out_dir = "./shared_volume/realtime_" + now +"/"
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
    obs_embeddings=["scvi"],
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
# ## Assign high-level cell types

# %%
print("assigning high-level cell types...")

def assign_high_level_cell_types(adata, high_level_cell_types):
    """
    Assign high-level cell types to cells in an AnnData object based on ontology ancestors.
    
    Parameters
    ----------
    adata : AnnData
        The AnnData object to annotate. Modified in place.
    high_level_cell_types : set or list
        A collection of high-level cell type names to match against ancestors.
    
    Returns
    -------
    None
        The function modifies adata.obs in place, adding two columns:
        - 'high_level_cell_type_ontology_term_id'
        - 'high_level_cell_type'
    """
    # Get all unique cell_type_ontology_term_id values except "unknown"
    cell_type_ids = adata.obs["cell_type_ontology_term_id"].unique()
    cell_type_ids = [ctid for ctid in cell_type_ids if ctid != "unknown"]
    
    # Map from cell_type_ontology_term_id to (high_level_id, high_level_name)
    adapter = get_adapter("sqlite:obo:cl")
    high_level_map = {}
    
    for ctid in cell_type_ids:
        try:
            ancestors = list(adapter.ancestors([ctid]))
        except Exception as e:
            print(f"Warning: Could not get ancestors for {ctid}: {e}")
            high_level_map[ctid] = (None, None)
            continue
        matching_ancestors = []
        for ancestor in ancestors:
            if str(ancestor).startswith("CL:"):
                name = adapter.label(ancestor)
                if name in high_level_cell_types:
                    matching_ancestors.append((ancestor, name))
        if len(matching_ancestors) == 1:
            ancestor, name = matching_ancestors[0]
            high_level_map[ctid] = (ancestor, name)
        elif len(matching_ancestors) == 0:
            print(f"Warning: No matching ancestor in high_level_cell_types for {ctid} ({adapter.label(ctid)})")
            high_level_map[ctid] = (None, None)
        elif (
            len(matching_ancestors) == 2 and
            ("hematopoietic cell" in [n for _, n in matching_ancestors]) and
            ("connective tissue cell" in [n for _, n in matching_ancestors])
        ):
            # Special case: show only hematopoietic cell
            for ancestor, name in matching_ancestors:
                if name == "hematopoietic cell":
                    high_level_map[ctid] = (ancestor, name)
                    break
        else:
            print(f"Warning: Multiple matching ancestors for {ctid} ({adapter.label(ctid)}): {matching_ancestors}")
            high_level_map[ctid] = (None, None)
    
    # Now, map these to the obs DataFrame
    def get_high_level_id(ctid):
        return high_level_map.get(ctid, (None, None))[0]
    
    def get_high_level_name(ctid):
        return high_level_map.get(ctid, (None, None))[1]
    
    adata.obs["high_level_cell_type_ontology_term_id"] = adata.obs["cell_type_ontology_term_id"].map(get_high_level_id)
    adata.obs["high_level_cell_type"] = adata.obs["cell_type_ontology_term_id"].map(get_high_level_name)

# Assign high-level cell types to adata_query (as per the comment in CELL INDEX 19)
assign_high_level_cell_types(adata_query, high_level_cell_types)

print(adata_query.obs[["cell_type","high_level_cell_type"]].value_counts())

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
    print("Package versions have been written to conda-requirements.txt")
else:
    print("Not in a conda environment. Running pip freeze")
    result = subprocess.run(['pip', 'freeze'], capture_output=True, text=True)
    packages = result.stdout
    # Write the packages to a file
    with open(out_dir + 'sys-requirements.txt', 'w') as f:
        f.write(packages)
os.chmod(out_dir, 0o777)
print("Done.")