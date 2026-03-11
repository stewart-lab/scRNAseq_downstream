# %% [markdown]
# # Annotate an scRNA-seq dataset using CASSIA

# %%
# ### Load libraries
# Standard python libraries
import numpy as np
import pandas as pd
import os

# scVerse
import anndata as ad
import scanpy as sc

# CASSIA
import CASSIA

# Local utilities
from utils import load_config, get_data_dir, initialize_output_directory, save_package_versions, compare_cell_metadata_cols, plot_umaps

# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
import warnings
from anndata import ImplicitModificationWarning
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

# %% [markdown]
# ### Specialized functions for this module
# <!-- Add any specialized functions for this module here. -->
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

    # Input files
    DATA_DIR = get_data_dir(config_dict)
    ref_data_file = DATA_DIR + config_dict[METHOD]["ref_data_file"]
    query_data_file = DATA_DIR + config_dict[METHOD]["query_data_file"]
    model_folder = DATA_DIR + config_dict[METHOD]["model_folder"]

    # Metadata column names
    gene_id_column = config_dict[METHOD]["gene_id_column"]
    cell_type_column = config_dict[METHOD]["cell_type_column"]
    high_level_cell_type_column = config_dict[METHOD]["high_level_cell_type_column"]

    # Output file
    out_dir = initialize_output_directory(config_dict)
    output_file = config_dict[METHOD]["output_file"]

    # Random seed for reproducibility
    random_seed = config_dict[METHOD].get("random_seed", 67)
    np.random.seed(random_seed)

    # Compare predicted cell types to gold standard annotations in the query dataset, if available
    compare_to_gold_standard = config_dict[METHOD].get("compare_to_gold_standard", False)
    gold_standard_cell_type_column = config_dict[METHOD].get("gold_standard_cell_type_column", None)
    gold_standard_high_level_cell_type_column = config_dict[METHOD].get("gold_standard_high_level_cell_type_column", None)

    # Annotation method-specific parameters
    # <!-- Add any additional method-specific parameters here. -->

    # ## Load data from files
    # Query dataset
    print()
    print(f"Loading query dataset from file: {query_data_file}")
    adata_query = sc.read_h5ad(query_data_file)
    print(adata_query)

    # Reference dataset
    print()
    print(f"Loading reference dataset from file: {ref_data_file}")
    adata_ref = sc.read_h5ad(ref_data_file)
    print(adata_ref)

    # Index genes by Ensembl ids
    adata_query.var.index = adata_query.var[gene_id_column]
    adata_ref.var.index = adata_ref.var[gene_id_column]
    print(adata_query.var)

    # ## Annotate the query dataset
    # <!-- Add your annotation code here. -->

    # ## Save the generated dataset and package versions

    # Save the annotated query dataset
    print()
    print("Saving the annotated query dataset...")
    adata_query.write_h5ad(os.path.join(out_dir, output_file))

    # Save package versions
    save_package_versions(out_dir)

    # Set permissions for the output directory to be readable and writable by all users
    os.chmod(out_dir, 0o777)

    # All done
    print()
    print("Done.")

# %% Run the main function
if __name__ == "__main__":
    main()