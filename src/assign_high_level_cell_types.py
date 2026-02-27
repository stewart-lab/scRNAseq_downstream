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
import os
import pandas as pd
import warnings

# scVerse
import anndata as ad
from anndata import ImplicitModificationWarning
import scanpy as sc

# Ontology
from oaklib import get_adapter

# Local utilities
from utils import load_config, get_data_dir, initialize_output_directory, save_package_versions

# Some settings to avoid errors/warnings
# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# %% [markdown]
# ## Specialized functions for this module

# %% Assign high-level cell types to cells in an AnnData object based on ontology ancestors.
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

# %% [markdown]
# # Main script
# %%
def main():
    # ## Parse the config file and set parameters
    # Configuration dictionary and input and output directories
    # Load configuration from config.json file
    config_dict = load_config()

    # Set input data directory
    DATA_DIR = get_data_dir(config_dict)

    # Initialie output directory
    out_dir = initialize_output_directory(config_dict)

    # Set custom parameters for the method
    # Method name
    METHOD = config_dict["METHOD"]

    # Input file
    input_file = DATA_DIR + config_dict[METHOD]["input_file"]

    # High-level cell types
    high_level_cell_types = (
        config_dict[METHOD]["high_level_cell_types"]
    )

    # Output file
    output_file = config_dict[METHOD]["output_file"]

    # ### Load the input dataset
    print(f"Loading input dataset from file: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(adata)

    # ## Assign high-level cell types
    print("assigning high-level cell types...")
    # Assign high-level cell types to adata
    assign_high_level_cell_types(adata, high_level_cell_types)

    # Print value counts for high and low level cell types
    counts = adata.obs[["high_level_cell_type","cell_type"]].value_counts()
    counts_sorted = (
        counts.rename("count")
        .reset_index()
        .query("count > 0")
        .sort_values(
            ["high_level_cell_type", "count", "cell_type"],
            ascending=[True, False, True],
        )
    )
    print(counts_sorted)

    # ## Save the generated dataset and package versions

    # Save the generated dataset
    print(f"Saving generated dataset to {output_file} ...")
    adata.write_h5ad(out_dir + output_file)

    # Save package versions
    save_package_versions(out_dir)

    # Set permissions for the output directory to be readable and writable by all users
    os.chmod(out_dir, 0o777)

    # All done
    print("Done.")

if __name__ == "__main__":
    main()