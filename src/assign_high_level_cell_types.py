# %% [markdown]
# # Assign high-level cell types based on existing low-level cell type annotations

# %% [markdown]
# This tool takes an annotated dataset and assigns high-level cell types based on a list supplied by the user. Existing low-level cell types, such as those included in CellxGene datasets, are mapped to a smaller set of high-level types using OBO Cell Ontology (CL). Both the existing cell type annotations and the high-level cell types must match valid ontology terms. The reason one may want to assign higher-level cell types is that they can be annotated with higher confidence.

# Pre-install the cellxgene_scvi conda environment from `cellxgene_scvi.yml` before attempting to run this tool.


# %% [markdown]
# ## Set up

# %% [markdown]
# ### Load libraries

# %%
# Standard python libraries
from datetime import datetime
import json
import os
import pandas as pd
import shutil
import subprocess
import warnings

# scVerse
import anndata as ad
from anndata import ImplicitModificationWarning
import scanpy as sc

# Ontology
from oaklib import get_adapter

# Some settings to avoid errors/warnings
# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# %% [markdown]
# ### Parse the config file.

# %%
# Read in the config file
work_dir = os.getcwd()
with open(work_dir+'/config.json') as f:
    config_dict = json.load(f)
    print(
        "loaded config file: ", 
        config_dict["assign_high_level_cell_types"]
    )

# docker and data directory
# (the data directory is not really used in this script, but we parse it for 
# consistency with other pipeline scripts)
docker = config_dict["docker"]
if docker == "TRUE" or docker == "true" or docker == "T" or docker == "t":
    DATA_DIR = "./data/input_data/"
else:
    DATA_DIR = config_dict["assign_high_level_cell_types"]["DATA_DIR"]

# Input file
input_file = config_dict["assign_high_level_cell_types"]["input_file"]

# High-level cell types
high_level_cell_types = (
    config_dict["assign_high_level_cell_types"]["high_level_cell_types"]
)

# Output file
output_file = config_dict["assign_high_level_cell_types"]["output_file"]

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
# ### Load the input dataset
print(f"Loading input dataset from file: {input_file}")
adata = sc.read_h5ad(input_file)
print(adata)

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

# Assign high-level cell types to adata and print the resulting value counts
assign_high_level_cell_types(adata, high_level_cell_types)
print(adata.obs[["cell_type","high_level_cell_type"]].value_counts())

# %% [markdown]
# ## Save the generated dataset and package versions

# %%
# Save the generated dataset
print(f"Saving generated dataset to {output_file} ...")
adata.write_h5ad(out_dir + output_file)

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