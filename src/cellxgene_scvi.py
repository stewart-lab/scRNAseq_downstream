# %% [markdown]
# # Annotate a dataset using scVI with a reference dataset from CellxGene

# %% [markdown]
# ## Set up

# %% [markdown]
# ### Load libraries

# %%
# Standard python libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import os

# scVerse
import anndata as ad
import scanpy as sc
import scvi

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

# Project a dataset to the scVI embedding
def project_adata_to_scvi(adata, model_folder):
    """
    Project query dataset into scVI embedding space using a pre-trained model.
    
    Parameters
    ----------
    adata : AnnData
        Query AnnData object to project.
    model_folder : str
        Path to the pre-trained scVI model folder. The folder must contain file "model.pt".
        
    Returns
    -------
    adata : AnnData
        Query AnnData object with scVI latent representation in obsm["scvi"].
    """
    # Copy the query dataset so as to preserve it. When the scVI model is 
    # loaded, all but 8,000 genes will be dropped.
    adata_scvi = adata.copy()

    # Load the scVI model and prepare the query data
    print()
    print("Preparing the query data for scVI model...")
    scvi.model.SCVI.prepare_query_anndata(adata_scvi, model_folder)
    print(adata_scvi)

    # Load the query data into the model, set "is_trained" to True to trick the 
    # model into thinking it was already trained, and do a forward pass through 
    # the model to get the latent representation of the query data.
    print()
    print("Projecting the query data to SCVI embedding...")
    vae_q = scvi.model.SCVI.load_query_data(
        adata_scvi,
        model_folder,
    )

    # This allows for a simple forward pass
    vae_q.is_trained = True
    latent = vae_q.get_latent_representation()
    adata.obsm["scvi"] = latent
    print(adata)
    
    return adata

# Predict cell types using a random forest classifier
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

# Compute fraction of correct predictions compared to gold standard annotations, excluding "unknown" annotations
def compute_frac_correct(
    adata,
    predicted_cell_type_column,
    gold_standard_cell_type_column,
):
    """
    Compute the fraction of correct predictions compared to gold standard annotations, excluding "unknown" annotations.
    Parameters
    ----------
    adata : AnnData
        AnnData object containing the predicted and gold standard cell type annotations in adata.obs.
    predicted_cell_type_column : str
        The name of the column in adata.obs containing the predicted cell type labels.
    gold_standard_cell_type_column : str
        The name of the column in adata.obs containing the gold standard cell type labels.

    Returns
    -------
    frac_correct : float
        The fraction of correct predictions compared to gold standard annotations, excluding "unknown" annotations.
    """
    predicted_cell_types = adata.obs[predicted_cell_type_column]
    gold_standard_cell_types = adata.obs[gold_standard_cell_type_column]
    mask = gold_standard_cell_types != "unknown"
    if mask.sum() == 0:
        return np.nan
    frac_correct = (
        predicted_cell_types[mask] == gold_standard_cell_types[mask]
    ).mean()
    return frac_correct

# Plot a set of UMAPs with consistent figure dimensions and spacing.
def plot_umaps(adata, out_dir, plot_specs, fixed_dpi=300, plot_area_in=5.2, left_margin_in=0.55, bottom_margin_in=0.55, top_margin_in=0.35, right_margin_continuous_in=2, right_margin_categorical_in=4):
    """
    Plot UMAPs with consistent figure dimensions and spacing.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object containing the data to plot.
    out_dir : str
        Directory to save the output plots.
    plot_specs : list of tuples
        List of (color_col, filename, title) tuples specifying what to plot.
    fixed_dpi : int, default 300
        DPI for saved figures.
    plot_area_in : float, default 5.2
        Width and height of the plot area in inches.
    left_margin_in : float, default 0.55
        Left margin in inches.
    bottom_margin_in : float, default 0.55
        Bottom margin in inches.
    top_margin_in : float, default 0.35
        Top margin in inches.
    right_margin_continuous_in : float, default 2
        Right margin in inches for continuous data legends.
    right_margin_categorical_in : float, default 4
        Right margin in inches for categorical legends.
    """

    for color_col, filename, title in plot_specs:
        series = adata.obs[color_col]
        is_continuous = pd.api.types.is_numeric_dtype(series)

        right_margin_in = (
            right_margin_continuous_in if is_continuous else right_margin_categorical_in
        )

        fig_w = left_margin_in + plot_area_in + right_margin_in
        fig_h = bottom_margin_in + plot_area_in + top_margin_in
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=fixed_dpi)

        # Position axis explicitly so the plotting area is always the same size
        ax_left = left_margin_in / fig_w
        ax_bottom = bottom_margin_in / fig_h
        ax_width = plot_area_in / fig_w
        ax_height = plot_area_in / fig_h

        with plt.rc_context():
            sc.pl.umap(
                adata,
                color=color_col,
                title=title,
                ax=ax,
                show=False,
                legend_loc="right margin",
                colorbar_loc="right",
            )

            ax.set_position([ax_left, ax_bottom, ax_width, ax_height])
            ax.set_box_aspect(1)

            fig.savefig(os.path.join(out_dir, filename), dpi=fixed_dpi, bbox_inches="tight", pad_inches=0.02)
            plt.close(fig)

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

    # ### Load data from files
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

    # ## Annotate the query dataset
    # Add necessary cell metadata columns
    print()
    print("Formatting the data...")
    adata_query.obs["n_counts"] = adata_query.X.sum(axis=1)
    adata_query.obs["joinid"] = list(range(adata_query.n_obs))
    adata_query.obs["batch"] = "unassigned"
    adata_query.obs["dataset_id"] = "query"

    adata_ref.obs["n_counts"] = adata_ref.X.sum(axis=1)
    adata_ref.obs["joinid"] = list(range(adata_ref.n_obs))
    adata_ref.obs["batch"] = "unassigned"
    adata_ref.obs["dataset_id"] = "reference"

    # Index genes by Ensembl ids
    adata_query.var.index = adata_query.var[gene_id_column]
    adata_ref.var.index = adata_ref.var[gene_id_column]
    print(adata_query.var)

    # ### Project the query and reference datasets into the scVI embedding
    # Project the query dataset into the scVI embedding
    adata_query = project_adata_to_scvi(adata_query, model_folder)
    adata_ref = project_adata_to_scvi(adata_ref, model_folder)

    # ### Predict high-level cell types
    # Predict high-level cell types in the query dataset
    print()
    print("Predicting high-level cell types for the query dataset...")
    preds, probs = predict_cell_types_with_rf(adata_query, adata_ref, high_level_cell_type_column)
    adata_query.obs["predicted_" + high_level_cell_type_column] = preds
    adata_query.obs[
        "predicted_" + high_level_cell_type_column + "_probability"
    ] = probs
    print()
    print("Probability distribution for predicted high-level cell types:")
    print(
        adata_query.obs[
            "predicted_" + high_level_cell_type_column + "_probability"
        ].describe()
    )

    # Print the counts of predicted high-level cell types in the query dataset
    print()
    print("Counts of predicted high-level cell types:")
    print(adata_query.obs["predicted_" + high_level_cell_type_column].value_counts())

    # Compare predicted high-level cell types to gold standard annotations, if available
    if compare_to_gold_standard and gold_standard_high_level_cell_type_column in adata_query.obs:
        print()
        print("Comparing predicted high-level cell types to gold standard annotations...")
        frac_correct = compute_frac_correct(
            adata_query,
            "predicted_" + high_level_cell_type_column,
            gold_standard_high_level_cell_type_column,
        )
        print(f"Fraction of correctly predicted high-level cell types (excluding 'unknown'): {frac_correct:.3f}")
        compare_cell_metadata_cols(
            "predicted_" + high_level_cell_type_column,
            gold_standard_high_level_cell_type_column,
            adata_query,
            out_dir
        )

    # ### Predict low-level cell types for each high-level cell type
    # For each high-level cell type, annotate within that group
    print()
    print("Predicting low-level cell types for each high-level cell type...")
    col_prev_pred = "predicted_" + high_level_cell_type_column
    unique_types = adata_query.obs[col_prev_pred].dropna().unique()
    for ct in unique_types:
        print()
        print(f"  Annotating {ct}...")
        # Mask for query and reference cells of this type at previous level
        mask_query = adata_query.obs[col_prev_pred] == ct
        mask_ref = adata_ref.obs[high_level_cell_type_column] == ct
        if mask_ref.sum() == 0 or mask_query.sum() == 0:
            print(f"    Skipping {ct}: no cells in reference or query.")
            continue
        else:
            print(f"    Reference cells: {mask_ref.sum()}, Query cells: {mask_query.sum()}")
        preds, probs = predict_cell_types_with_rf(
            adata_query, adata_ref, cell_type_column, mask_query=mask_query, mask_ref=mask_ref
        )
        # Assign only to the relevant subset
        adata_query.obs.loc[mask_query, "predicted_" + cell_type_column] = preds
        adata_query.obs.loc[mask_query, "predicted_" + cell_type_column + "_probability"] = probs
    
    # Print counts of predicted cell types
    print()
    print("Counts of predicted cell types:")
    print(
        adata_query.obs[[
            "predicted_" + high_level_cell_type_column,
            "predicted_" + cell_type_column
        ]].value_counts()
    )

    # Print the distribution of prediction probabilities for the low-level cell types
    print()
    print("Probability distribution for predicted low-level cell types:")
    print(adata_query.obs["predicted_" + cell_type_column + "_probability"].describe())

    # Compare predicted low-level cell types to gold standard annotations, if available
    if compare_to_gold_standard and gold_standard_cell_type_column in adata_query.obs:
        print()
        print("Comparing predicted low-level cell types to gold standard annotations...")
        frac_correct = compute_frac_correct(
            adata_query,
            "predicted_" + cell_type_column,
            gold_standard_cell_type_column,
        )
        print(f"Fraction of correctly predicted low-level cell types (excluding 'unknown'): {frac_correct:.3f}")
        compare_cell_metadata_cols(
            "predicted_" + cell_type_column,
            gold_standard_cell_type_column,
            adata_query,
            out_dir
        )

    # ## Visualize the results
    # ### High-level cell types
    print()
    print("Combining query and reference datasets for visualization...")
    # Set "predicted" cell types in the reference dataset to the actual reference cell types
    adata_ref.obs["predicted_" + high_level_cell_type_column] = adata_ref.obs[high_level_cell_type_column]
    adata_ref.obs["predicted_" + cell_type_column] = adata_ref.obs[cell_type_column]

    # Set the reference prediction scores to nan
    adata_ref.obs["predicted_" + high_level_cell_type_column + "_probability"] = np.nan
    adata_ref.obs["predicted_" + cell_type_column + "_probability"] = np.nan

    # Combine the datasets
    adata_combined = ad.concat([adata_query, adata_ref])

    # Plot the UMAPs
    print("Plotting UMAPs with predicted high-level cell types...")
    sc.pp.neighbors(adata_combined, n_neighbors=15, use_rep="scvi", metric="correlation")
    sc.tl.umap(adata_combined)
    plot_specs = [
        (
            "dataset_id", 
            "dataset_id_umap.pdf", 
            "Dataset"
        ),
        (
            "predicted_" + high_level_cell_type_column,
            "predicted_high_level_cell_type_umap.pdf",
            "Predicted High-Level Cell Type",
        ),
        (
            "predicted_" + high_level_cell_type_column + "_probability",
            "predicted_high_level_cell_type_probability_umap.pdf",
            "Prediction Probability",
        ),
    ]
    plot_umaps(adata_combined, out_dir, plot_specs)

    # ### Low-level cell types
    print()
    print("Plotting UMAPs for low-level cell types within each high-level cell type...")
    high_level_cell_types = adata_combined.obs["predicted_" + high_level_cell_type_column].unique()
    for hlct in high_level_cell_types:
        # Print the low-level cell types present in this high-level cell type
        print(f"### High-level cell type: {hlct}")
        mask = adata_combined.obs["predicted_" + high_level_cell_type_column] == hlct
        cell_types = adata_combined.obs.loc[mask, "predicted_" + cell_type_column].unique()
        cell_types_nonempty = [ct for ct in cell_types if pd.notna(ct)]
        print(f"  Low-level cell types: {cell_types_nonempty}")

        # Create a "safe" version of the high-level cell type name for use in file names.
        safe_hlct = str(hlct).replace("/", "_").replace(" ", "_")

        # Subset the data to just this high-level cell type for plotting.
        subset = adata_combined[mask, :]

        # Plot the UMAPs colored by dataset, low-level cell type and prediction probability
        low_level_plot_specs = [
            (
                "dataset_id", 
                f"{safe_hlct}_dataset_umap.pdf", 
                f"{hlct} - Dataset"
            ),
            (
                "predicted_" + cell_type_column,
                f"{safe_hlct}_predicted_cell_type_umap.pdf",
                f"{hlct} - Predicted Cell Type",
            ),
            (
                "predicted_" + cell_type_column + "_probability",
                f"{safe_hlct}_prediction_probability_umap.pdf",
                f"{hlct} - Prediction Probability",
            ),
        ]
        plot_umaps(subset, out_dir, low_level_plot_specs)

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