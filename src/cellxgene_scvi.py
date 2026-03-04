# %% [markdown]
# # Annotate a dataset using scVI with a reference dataset from CellxGene

# %% [markdown]
# ## Set up

# %% [markdown]
# ### Load libraries

# %%
# Standard python libraries
from datetime import datetime
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shutil
from sklearn.ensemble import RandomForestClassifier
import subprocess
import os

# scVerse
import anndata as ad
import scanpy as sc
import scvi

# Local utilities
from utils import load_config, get_data_dir, initialize_output_directory, save_package_versions

# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
import warnings
from anndata import ImplicitModificationWarning
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

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
    model_file = DATA_DIR + config_dict[METHOD]["model_file"]

    # Output file
    out_dir = initialize_output_directory(config_dict)
    output_file = config_dict[METHOD]["output_file"]

    # Gene id column
    gene_id_col = config_dict[METHOD]["gene_id_col"]
    # CONTINUE HERE
    # add gene_id_col to the config file and use it to set the var index for both the query and reference datasets
    # modify ids in the query dataset (use query dataset without added noise)

    # ### Load data from files
    # Query dataset
    print(f"Loading query dataset from file: {query_data_file}")
    adata_query = sc.read_h5ad(query_data_file)
    print(adata_query)

    # Reference dataset
    print(f"Loading reference dataset from file: {ref_data_file}")
    adata_ref = sc.read_h5ad(ref_data_file)
    print(adata_ref)

    # ## Annotate the query dataset
    # Add necessary cell metadata columns
    print("Formatting the data...")
    adata_query.obs["n_counts"] = adata_query.X.sum(axis=1)
    adata_query.obs["joinid"] = list(range(adata_query.n_obs))
    adata_query.obs["batch"] = "unassigned"
    adata_query.obs["dataset_id"] = "query"

    adata_ref.obs["n_counts"] = adata_ref.X.sum(axis=1)
    adata_ref.obs["joinid"] = list(range(adata_ref.n_obs))
    adata_ref.obs["batch"] = "unassigned"
    adata_ref.obs["dataset_id"] = "reference"

    # %% [markdown]
    # Index genes by Ensembl ids

    # %%
    adata_query.var.index = adata_query.var["feature_id"]
    adata_ref.var.index = adata_ref.var["feature_id"]
    print(adata_query.var)

    # %% [markdown]
    # ### Project the query dataset into the scVI embedding

    # %% [markdown]
    # Copy the query dataset so as to preserve it. When the scVI model is loaded, all but 8,000 genes will be dropped.

    # %%
    adata_query_scvi = adata_query.copy()

    # %% [markdown]
    # Load the scVI model and prepare the query data

    # %%
    print("Preparing the query data for scVI model...")
    scvi.model.SCVI.prepare_query_anndata(adata_query_scvi, model_file)
    print(adata_query_scvi)

    # %% [markdown]
    # Load the query data into the model, set “is_trained” to True to trick the model into thinking it was already trained, and do a forward pass through the model to get the latent representation of the query data.

    # %%
    print("Projecting the query data to SCVI embedding...")
    vae_q = scvi.model.SCVI.load_query_data(
        adata_query_scvi,
        model_file,
    )

    # This allows for a simple forward pass
    vae_q.is_trained = True
    latent = vae_q.get_latent_representation()
    adata_query.obsm["scvi"] = latent
    print(adata_query)

    # %% [markdown]
    # Combine and plot the two datasets

    # %%
    print("Combining the query and reference datasets for visualization...")
    adata_combined = ad.concat([adata_query, adata_ref])
    adata_combined.obs_names_make_unique() # this is necessary when the query and the reference have been sampled from the same parent dataset
    sc.pp.neighbors(
        adata_combined, n_neighbors=15, use_rep="scvi", metric="correlation"
    )
    sc.tl.umap(adata_combined)
    with plt.rc_context():
        sc.pl.umap(adata_combined, color=["dataset_id"])
        plt.savefig(out_dir + "query_and_reference_umap.pdf")

    # %% [markdown]
    # ### Predict high-level cell types

    # %%
    def predict_cell_types_with_rf(
        adata_query,
        adata_ref,
        cell_type_col,
        mask_query=None,
        mask_ref=None
    ):
        """
        Fit a Random Forest Classifier on the reference dataset's embedding and use it to predict cell type labels
        for the query dataset, for the specified cell type column.
        Also computes and plots the prediction confidence probabilities.

        Parameters
        ----------
        adata_query : AnnData
            Query AnnData object to annotate.
        adata_ref : AnnData
            Reference AnnData object with cell type annotations.
        cell_type_col : str
            The name of the cell type column to use for training and prediction.
        mask_query : pd.Series, np.ndarray, or None
            Boolean mask for adata_query.obs to select cells to annotate. If None, use all.
        mask_ref : pd.Series, np.ndarray, or None
            Boolean mask for adata_ref.obs to select reference cells for training. If None, use all.

        Returns
        -------
        preds : np.ndarray
            Predicted cell type labels for the selected query cells.
        probs : np.ndarray
            Prediction confidence probabilities for the selected query cells.
        """

        # Select reference cells
        if mask_ref is not None:
            X_ref = adata_ref.obsm["scvi"][mask_ref]
            y_ref = adata_ref.obs[cell_type_col][mask_ref].values
        else:
            X_ref = adata_ref.obsm["scvi"]
            y_ref = adata_ref.obs[cell_type_col].values

        # Select query cells
        if mask_query is not None:
            X_query = adata_query.obsm["scvi"][mask_query]
        else:
            X_query = adata_query.obsm["scvi"]

        # Fit classifier
        rfc = RandomForestClassifier()
        rfc.fit(X_ref, y_ref)
        preds = rfc.predict(X_query)

        # Compute confidence scores
        probabilities = rfc.predict_proba(X_query)
        confidence = np.zeros(len(preds))
        for i in range(len(preds)):
            confidence[i] = probabilities[i][rfc.classes_ == preds[i]].item()

        return preds, confidence

    # %% [markdown]
    # Predict high-level cell types in the query dataset

    # %%
    print("Predicting high-level cell types for the query dataset...")
    preds, probs = predict_cell_types_with_rf(adata_query, adata_ref, "high_level_cell_type")
    adata_query.obs["predicted_high_level_cell_type"] = preds
    adata_query.obs["predicted_high_level_cell_type_probability"] = probs
    print("Probability distribution for predicted high-level cell types:")
    print(adata_query.obs["predicted_high_level_cell_type_probability"].describe())

    # %%
    print("Counts of predicted high-level cell types:")
    print(adata_query.obs["predicted_high_level_cell_type"].value_counts())

    # %% [markdown]
    # ### Predict low-level cell types for each high-level cell type

    # %%
    # For each high-level cell type, annotate within that group
    print("Predicting low-level cell types for each high-level cell type...")
    col_prev_pred = "predicted_high_level_cell_type"
    unique_types = adata_query.obs[col_prev_pred].dropna().unique()
    for ct in unique_types:
        print(f"  Annotating {ct}...")
        # Mask for query and reference cells of this type at previous level
        mask_query = adata_query.obs[col_prev_pred] == ct
        mask_ref = adata_ref.obs["high_level_cell_type"] == ct
        if mask_ref.sum() == 0 or mask_query.sum() == 0:
            print(f"    Skipping {ct}: no cells in reference or query.")
            continue
        else:
            print(f"    Reference cells: {mask_ref.sum()}, Query cells: {mask_query.sum()}")
        preds, probs = predict_cell_types_with_rf(
            adata_query, adata_ref, "cell_type", mask_query=mask_query, mask_ref=mask_ref
        )
        # Assign only to the relevant subset
        adata_query.obs.loc[mask_query, "predicted_cell_type"] = preds
        adata_query.obs.loc[mask_query, "predicted_cell_type_probability"] = probs

    # %%
    print("Counts of predicted cell types:")
    print(
        adata_query.obs[[
            "predicted_high_level_cell_type",
            "predicted_cell_type"
        ]].value_counts()
    )

    # %%
    print(adata_query.obs["predicted_cell_type_probability"].describe())

    # %% [markdown]
    # ## Visualize the results

    # %% [markdown]
    # ### High-level cell types

    # %%
    print("Combining query and reference datasets for visualization...")
    # Set "predicted" cell types in the reference dataset to the actual reference cell types
    adata_ref.obs["predicted_high_level_cell_type"] = adata_ref.obs["high_level_cell_type"]
    adata_ref.obs["predicted_cell_type"] = adata_ref.obs["cell_type"]

    # Set the reference prediction scores to nan
    adata_ref.obs["predicted_high_level_cell_type_probability"] = np.nan
    adata_ref.obs["predicted_cell_type_probability"] = np.nan

    # Combine the datasets
    adata_combined = ad.concat([adata_query, adata_ref])
    adata_combined.obs_names_make_unique() # this is necessary when the query and the reference have been sampled from the same parent dataset

    # Plot the UMAP
    print("Plotting UMAP with predicted high-level cell types...")
    sc.pp.neighbors(adata_combined, n_neighbors=15, use_rep="scvi", metric="correlation")
    sc.tl.umap(adata_combined)
    with plt.rc_context():
        sc.pl.umap(
            adata_combined, 
            color=[
                "dataset_id", 
                "predicted_high_level_cell_type",
                "predicted_high_level_cell_type_probability"
            ]
        )
        plt.savefig(out_dir + "high_level_cell_types_umap.pdf")

    # %% [markdown]
    # ### Low-level cell types

    # %%
    print("Plotting UMAPs for low-level cell types within each high-level cell type...")
    high_level_cell_types = adata_combined.obs["predicted_high_level_cell_type"].unique()
    for hlct in high_level_cell_types:
        print(f"### High-level cell type: {hlct}")
        mask = adata_combined.obs["predicted_high_level_cell_type"] == hlct
        cell_types = adata_combined.obs.loc[mask, "predicted_cell_type"].unique()
        cell_types_nonempty = [ct for ct in cell_types if pd.notna(ct)]
        print(f"  Low-level cell types: {cell_types_nonempty}")
        with plt.rc_context():
            sc.pl.umap(
                adata_combined[mask, :],
                color=[
                    "dataset_id", 
                    "predicted_cell_type",
                    "predicted_cell_type_probability"
                ],
                title=[
                    f"{hlct} - Dataset", 
                    f"{hlct} - Predicted Cell Type", 
                    f"{hlct} - Prediction Probability"
                ]
            )
            plt.savefig(out_dir + f"{hlct}_low_level_cell_types_umap.pdf")   

    # %% [markdown]
    # ## Save the annotated query dataset

    # %%
    print("Saving the annotated query dataset...")
    adata_query.write_h5ad(out_dir + output_file)

    #%%
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