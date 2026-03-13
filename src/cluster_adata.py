# %% [markdown]
# # Annotate an scRNA-seq dataset using ...

# %%
# ### Load libraries
# Standard python libraries
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

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
# Compute and save DE genes
def compute_and_save_cluster_de(
    adata: ad.AnnData,
    cluster_column: str,
    de_method: str,
    key_added: str,
    gene_symbol_column: str,
    out_dir: str,
    output_file_dotplot: str,
    output_file_de: str,
):
    """
    Compute and save cluster-specific differentially expressed genes.
    
    Parameters    ----------
    adata: anndata.AnnData
        The annotated dataset containing the clusters for which to compute DE genes.
    cluster_column: str
        The name of the column in adata.obs containing the cluster labels.
    de_method: str
        The method to use for computing DE genes (e.g., "wilcoxon", "t-test", etc.).
    key_added: str
        The key to add to adata.uns to store the DE results.
    gene_symbol_column: str
        The name of the column in adata.var containing gene symbols.
    out_dir: str
        The directory in which to save the DE results.
    output_file_dotplot: str
        The name of the file in which to save the DE dotplot.
    output_file_de: str
        The name of the file in which to save the DE results as a CSV file.

    Returns    -------
    None
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_column,
        method=de_method,
        key_added=key_added,
    )

    sc.tl.dendrogram(adata, groupby=cluster_column, use_rep="scvi")

    plot_gene_symbol_column = gene_symbol_column
    if gene_symbol_column in adata.var.columns:
        duplicate_symbol_mask = adata.var[gene_symbol_column].astype("string").duplicated(keep=False)
        if duplicate_symbol_mask.any():
            plot_gene_symbol_column = f"{gene_symbol_column}_unique_for_plot"
            adata.var[plot_gene_symbol_column] = adata.var[gene_symbol_column].astype("string")
            adata.var.loc[duplicate_symbol_mask, plot_gene_symbol_column] = (
                adata.var.loc[duplicate_symbol_mask, gene_symbol_column].astype("string")
                + " ("
                + adata.var.loc[duplicate_symbol_mask].index.astype("string")
                + ")"
            )

    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=cluster_column,
        standard_scale="var",
        n_genes=5,
        gene_symbols=plot_gene_symbol_column,
        key=key_added,
        show=False,
    )
    plt.gcf().savefig(os.path.join(out_dir, output_file_dotplot))
    plt.close()

    de_clusters_df = sc.get.rank_genes_groups_df(
        adata,
        group=None,
        key=key_added,
    )
    de_clusters_df.to_csv(os.path.join(out_dir, output_file_de), index=False)
    print(de_clusters_df)

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

    # Gene ids (Ensembl) and symbols
    gene_id_column = config_dict[METHOD]["gene_id_column"]
    gene_symbol_column = config_dict[METHOD]["gene_symbol_column"]

    # Output file
    out_dir = initialize_output_directory(config_dict)
    output_file_prefix = config_dict[METHOD]["output_file_prefix"]
    
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

    # Differential expression parameters
    de_method = config_dict[METHOD].get("de_method", "wilcoxon")

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
    sc.tl.leiden(adata, resolution=clustering_resolution, key_added=cluster_column, flavor=leiden_flavor, n_iterations=num_iterations, random_state=random_seed)
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
        cluster_mask = (
            (adata.obs[cluster_column].astype("string") == cluster_id)
            .fillna(False)
            .to_numpy(dtype=bool)
        )
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
            random_state=random_seed,
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
            compare_cell_metadata_cols(cluster_column, gold_standard_high_level_cell_type_column, adata, out_dir)

        if gold_standard_cell_type_column is not None:
            print()
            print(f"Comparing sub-clusters to gold standard cell type annotations (column: {gold_standard_cell_type_column})...")
            compare_cell_metadata_cols(subcluster_column, gold_standard_cell_type_column, adata, out_dir)

    # ## Compute differential expression between clusters
    # Normalize the data before computing differential expression
    print()
    print("Normalizing the data before computing differential expression...")

    # Save count data
    adata.layers["counts"] = adata.X.copy()

    # Normalize to median total counts
    sc.pp.normalize_total(adata)

    # Logarithmize the data
    sc.pp.log1p(adata)
    adata

    # Compute cluster-specific differentially expressed genes
    print()
    print("Computing differential expression between clusters...")
    compute_and_save_cluster_de(
        adata=adata,
        cluster_column=cluster_column,
        de_method=de_method,
        key_added="de_clusters",
        gene_symbol_column=gene_symbol_column,
        out_dir=out_dir,
        output_file_dotplot=output_file_prefix + "_cluster_de_dotplot.png",
        output_file_de=output_file_prefix + "_cluster_de.csv",
    )

    for cluster_id in cluster_ids:
        cluster_mask = (
            (adata.obs[cluster_column].astype("string") == cluster_id)
            .fillna(False)
            .to_numpy(dtype=bool)
        )
        adata_cluster = adata[cluster_mask].copy()

        n_subclusters = (
            adata_cluster.obs[subcluster_column]
            .astype("string")
            .dropna()
            .nunique()
        )

        if n_subclusters < 2:
            continue

        safe_cluster_id = str(cluster_id).replace(os.sep, "_")

        compute_and_save_cluster_de(
            adata=adata_cluster,
            cluster_column=subcluster_column,
            de_method=de_method,
            key_added=f"de_subclusters_{safe_cluster_id}",
            gene_symbol_column=gene_symbol_column,
            out_dir=out_dir,
            output_file_dotplot=f"{output_file_prefix}_{safe_cluster_id}_subcluster_de_dotplot.png",
            output_file_de=f"{output_file_prefix}_{safe_cluster_id}_subcluster_de.csv",
        )

    # ## Save the generated dataset and package versions

    # Save the clustered dataset
    print()
    print("Saving the clustered dataset...")
    adata.write_h5ad(os.path.join(out_dir, output_file_prefix + ".h5ad"))

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