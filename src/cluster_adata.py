# %% [markdown]
# # Annotate an scRNA-seq dataset using ...

# %%
# ### Load libraries
# Standard python libraries
from sklearn.metrics._scorer import metric
import numpy as np
import pandas as pd
import os

# scVerse
import anndata as ad
import scanpy as sc

# Local utilities
from utils import load_config, get_data_dir, initialize_output_directory, save_package_versions, compare_cell_metadata_cols

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

    # Input file
    DATA_DIR = get_data_dir(config_dict)
    input_data_file = DATA_DIR + config_dict[METHOD]["input_data_file"]

    # Gene ids (Ensembl)
    gene_id_column = config_dict[METHOD]["gene_id_column"]

    # Output file
    out_dir = initialize_output_directory(config_dict)
    output_file = config_dict[METHOD]["output_file"]

    # Random seed for reproducibility
    random_seed = config_dict[METHOD].get("random_seed", 67)
    np.random.seed(random_seed)

    # Compare clusters to gold standard annotations, if available
    compare_to_gold_standard = config_dict[METHOD].get("compare_to_gold_standard", False)
    gold_standard_high_level_cell_type_column = config_dict[METHOD].get("gold_standard_high_level_cell_type_column", None)
    gold_standard_cell_type_column = config_dict[METHOD].get("gold_standard_cell_type_column", None)

    # ### Annotation method-specific parameters
    # Neighborhood graph parameters
    embedding = config_dict[METHOD].get("embedding", "scvi")
    n_neighbors = config_dict[METHOD].get("n_neighbors", 15)
    distance_metric = config_dict[METHOD].get("distance_metric", "correlation")

    # Clustering parameters
    clustering_resolution = config_dict[METHOD].get("clustering_resolution", 0.1)
    subclustering_resolution = config_dict[METHOD].get("subclustering_resolution", 1.0)
    leiden_flavor = config_dict[METHOD].get("leiden_flavor", "igraph")
    num_iterations = config_dict[METHOD].get("num_iterations", 2)

    # Metadata column names for computed clusters and subclusters
    cluster_column = config_dict[METHOD]["cluster_column"]
    subcluster_column = config_dict[METHOD]["subcluster_column"]

    # ## Load data from file
    print()
    print(f"Loading input dataset from file: {input_data_file}")
    adata = sc.read_h5ad(input_data_file)
    print(adata)

    # Index genes by Ensembl ids
    adata.var.index = adata.var[gene_id_column]
    print()
    print("Gene metadata after indexing by gene ids:")
    print(adata.var)

    # ## Cluster the dataset
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=embedding, metric=distance_metric)

    # Cluster the dataset using Leiden clustering
    sc.tl.leiden(adata, resolution=clustering_resolution, key_added=cluster_column, flavor=leiden_flavor)
    print()
    print(f"Clustering results (column: {cluster_column}):")
    print(adata.obs[cluster_column].value_counts())

    # Subcluster each cluster
    adata.obs[subcluster_column] = pd.Series(index=adata.obs_names, dtype="string")

    cluster_ids = (
        adata.obs[cluster_column]
        .astype("string")
        .dropna()
        .unique()
        .tolist()
    )

    for cluster_id in cluster_ids:
        cluster_mask = adata.obs[cluster_column].astype("string") == cluster_id
        adata_cluster = adata[cluster_mask].copy()

        if adata_cluster.n_obs < 2:
            adata.obs.loc[adata_cluster.obs_names, subcluster_column] = f"{cluster_id}_0"
            continue

        cluster_n_neighbors = min(n_neighbors, adata_cluster.n_obs - 1)

        sc.pp.neighbors(
            adata_cluster,
            n_neighbors=cluster_n_neighbors,
            use_rep=embedding,
            metric=distance_metric,
        )

        sc.tl.leiden(
            adata_cluster,
            resolution=subclustering_resolution,
            key_added="_subcluster",
            flavor=leiden_flavor,
            n_iterations=num_iterations,
        )

        adata.obs.loc[adata_cluster.obs_names, subcluster_column] = (
            cluster_id + "_" + adata_cluster.obs["_subcluster"].astype("string")
        ).to_numpy()

    adata.obs[subcluster_column] = adata.obs[subcluster_column].astype("category")

    # ## Compare clusters to gold standard annotations, if available
    if compare_to_gold_standard:
        print()
        print("Comparing clusters to gold standard annotations...")

        if gold_standard_high_level_cell_type_column is not None:
            print()
            print(f"Comparing clusters to gold standard high-level cell type annotations (column: {gold_standard_high_level_cell_type_column})...")
            compare_cell_metadata_cols(adata, cluster_column, gold_standard_high_level_cell_type_column)

        if gold_standard_cell_type_column is not None:
            print()
            print(f"Comparing sub-clusters to gold standard cell type annotations (column: {gold_standard_cell_type_column})...")
            compare_cell_metadata_cols(adata, subcluster_column, gold_standard_cell_type_column)

    # ## Save the generated dataset and package versions

    # Save the clustered dataset
    print()
    print("Saving the clustered dataset...")
    adata.write_h5ad(os.path.join(out_dir, output_file))

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