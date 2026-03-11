# %% [markdown]
# # Test CASSIA on the Yayon dataset

# %% [markdown]
# ## Set up

# %% [markdown]
# ### Load libraries

# %%
import anndata
import cellxgene_census
import cellxgene_census.experimental
import numpy as np
import scanpy as sc
import scvi
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from oaklib import get_adapter
import CASSIA
import os

# %% [markdown]
# ### Load the Yayon dataset

# %% [markdown]
# Create a census object

# %%
census_version = "2025-11-08"
organism = "homo_sapiens"
census = cellxgene_census.open_soma(census_version=census_version)
census

# %% [markdown]
# Load the Yayon et al. dataset

# %%
dataset_ids = [
    "29244e1d-02e6-4133-b411-516ef7474638",
]

adata_census = cellxgene_census.get_anndata(
    census=census,
    measurement_name="RNA",
    organism="Homo sapiens",
    obs_value_filter=f"dataset_id in {dataset_ids}",
    obs_embeddings=["scvi"],
)
adata_census

# %% [markdown]
# ## Assign high-level cell types

# %% [markdown]
# Cell types and counts

# %%
adata_census.obs[["cell_type_ontology_term_id","cell_type"]].value_counts()

# %% [markdown]
# Drop the "unknown" cell type

# %%
# Steps:
# 1) record starting count
# 2) build mask excluding "unknown"
# 3) subset and copy AnnData
# 4) return dropped and remaining counts

before = int(adata_census.n_obs)
mask = adata_census.obs["cell_type_ontology_term_id"].fillna("") != "unknown"
adata_census = adata_census[mask].copy()
dropped = before - int(adata_census.n_obs)
{"dropped": dropped, "remaining": int(adata_census.n_obs)}

# %% [markdown]
# Connect to the Cell Ontology (CL)

# %%
adapter = get_adapter("sqlite:obo:cl")
adapter

# %% [markdown]
# Assign high-level cell types to all cells

# %%
# For each distinct cell_type_ontology_term_id (except "unknown"), find the ancestor from among target_labels,
# and save the ancestor's id and name in new adata_census.obs columns.

target_labels = {
    "epithelial cell", "hematopoietic cell", "neural cell", "connective tissue cell", "muscle cell",
    "absorptive cell", "endothelial cell", "transit amplifying cell", "embryonic stem cell", "pluripotent stem cell"
}

# Get all unique cell_type_ontology_term_id values except "unknown"
cell_type_ids = adata_census.obs["cell_type_ontology_term_id"].unique()
cell_type_ids = [ctid for ctid in cell_type_ids if ctid != "unknown"]

# Map from cell_type_ontology_term_id to (high_level_id, high_level_name)
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
            if name in target_labels:
                matching_ancestors.append((ancestor, name))
    if len(matching_ancestors) == 1:
        ancestor, name = matching_ancestors[0]
        high_level_map[ctid] = (ancestor, name)
    elif len(matching_ancestors) == 0:
        print(f"Warning: No matching ancestor in target_labels for {ctid} ({adapter.label(ctid)})")
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

adata_census.obs["high_level_cell_type_ontology_term_id"] = adata_census.obs["cell_type_ontology_term_id"].map(get_high_level_id)
adata_census.obs["high_level_cell_type"] = adata_census.obs["cell_type_ontology_term_id"].map(get_high_level_name)

adata_census.obs[["cell_type","high_level_cell_type"]].value_counts()

# %% [markdown]
# ## Prepare query and reference datasets

# %% [markdown]
# Create two random subsets having 100 representatives of each cell type

# %%
# Get all unique cell types
cell_types = adata_census.obs.loc[adata_census.obs['cell_type'] != "unknown", 'cell_type'].dropna().unique()

# For reproducibility
np.random.seed(86)

# For each cell type, sample up to 100 cells (if available)
indices_q = []
indices_ref = []

for ct in cell_types:
    idx = adata_census.obs[adata_census.obs['cell_type'] == ct].index
    n = min(100, len(idx))
    if n == 0:
        continue
    # Shuffle indices
    idx = np.random.permutation(idx)
    # If more than 200, split into two non-overlapping sets of 100
    if len(idx) >= 200:
        indices_q.extend(idx[:100])
        indices_ref.extend(idx[100:200])
    else:
        # If less than 200, split as evenly as possible
        split = n // 2
        indices_q.extend(idx[:split])
        indices_ref.extend(idx[split:n])
        
# If some cell types have <100 cells, the subsets will have fewer for those types

# Create the two AnnData subsets
adata_query = adata_census[indices_q, :].copy()
adata_ref = adata_census[indices_ref, :].copy()
adata_query, adata_ref

# %% [markdown]
# Tabulate the number of cells of each cell type in adata_query and adata_ref

# %%
# Tabulate the number of cells of each level 4 subtype in adata_query and adata_ref
tab_q = adata_query.obs['cell_type'].value_counts().sort_index()
tab_ref = adata_ref.obs['cell_type'].value_counts().sort_index()

# Combine into a DataFrame for easy comparison
celltype_counts = pd.DataFrame({
    "adata_query": tab_q,
    "adata_ref": tab_ref
}).fillna(0).astype(int)

celltype_counts

# %% [markdown]
# Add necessary cell metadata columns

# %%
adata_query.obs["n_counts"] = adata_query.X.sum(axis=1)
adata_query.obs["joinid"] = list(range(adata_query.n_obs))
adata_query.obs["batch"] = "unassigned"
adata_query.obs["dataset_id"] = "query"

adata_ref.obs["n_counts"] = adata_ref.X.sum(axis=1)
adata_ref.obs["joinid"] = list(range(adata_ref.n_obs))
adata_ref.obs["batch"] = "unassigned"
adata_ref.obs["dataset_id"] = "reference"

# %% [markdown]
# Combine and plot the two subsets

# %%
adata_combined = anndata.concat([adata_query, adata_ref])
sc.pp.neighbors(adata_combined, n_neighbors=15, use_rep="scvi", metric="correlation")
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=["dataset_id"])

# %% [markdown]
# Index genes by Ensembl ids

# %%
adata_query.var.index = adata_query.var["feature_id"]
adata_ref.var.index = adata_ref.var["feature_id"]
adata_query.var

# %% [markdown]
# ## Annotate the test dataset with scVI

# %% [markdown]
# ### Project the query dataset into the scVI embedding

# %% [markdown]
# Copy the query dataset so as to preserve it. When the scVI model is loaded, all but 8,000 genes will be dropped.

# %%
adata_query_scvi = adata_query.copy()
adata_query_scvi

# %% [markdown]
# Load the scVI model and prepare the query data (the model was downloaded from s3://cellxgene-contrib-public/models/scvi/2025-11-08/homo_sapiens/model.pt)

# %%
model_filename = "/mnt/cephfs/mir/rstewart/BrownCollab/BrownCollabData/CellxGene-scVI-models/2025-11-08-scvi-homo-sapiens"
scvi.model.SCVI.prepare_query_anndata(adata_query_scvi, model_filename)
adata_query_scvi

# %% [markdown]
# Load the query data into the model, set “is_trained” to True to trick the model into thinking it was already trained, and do a forward pass through the model to get the latent reprsentation of the query data.

# %%
vae_q = scvi.model.SCVI.load_query_data(
    adata_query_scvi,
    model_filename,
)

# This allows for a simple forward pass
vae_q.is_trained = True
latent = vae_q.get_latent_representation()
adata_query_scvi.obsm["scvi2"] = latent
adata_query_scvi

# %% [markdown]
# filter out missing features (not really necessary here)

# %%
adata_query_scvi = adata_query_scvi[:, adata_query_scvi.var["feature_name"].notnull().values].copy()
adata_query_scvi

# %% [markdown]
# Verify that the newly computed embedding is like the old one

# %%
adata_ref.obsm["scvi2"] = adata_ref.obsm["scvi"]
adata_combined = anndata.concat([adata_query_scvi, adata_ref])
sc.pp.neighbors(adata_combined, n_neighbors=15, use_rep="scvi2", metric="correlation")
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=["dataset_id"])

# %% [markdown]
# ### Predict high-level cell types

# %% [markdown]
# Tabulate high-level cell types in the reference dataset

# %%
adata_ref.obs["high_level_cell_type"].value_counts()

# %% [markdown]
# Define a function to annotate cell types based on a column in the reference

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
        X_ref = adata_ref.obsm["scvi2"][mask_ref]
        y_ref = adata_ref.obs[cell_type_col][mask_ref].values
    else:
        X_ref = adata_ref.obsm["scvi2"]
        y_ref = adata_ref.obs[cell_type_col].values

    # Select query cells
    if mask_query is not None:
        X_query = adata_query.obsm["scvi2"][mask_query]
    else:
        X_query = adata_query.obsm["scvi2"]

    # Fit classifier
    rfc = RandomForestClassifier()
    rfc.fit(X_ref, y_ref)
    preds = rfc.predict(X_query)

    # Compute confidence scores
    probabilities = rfc.predict_proba(X_query)
    confidence = np.zeros(len(preds))
    for i in range(len(preds)):
        confidence[i] = probabilities[i][rfc.classes_ == preds[i]]

    return preds, confidence

# %% [markdown]
# Predict high-level cell types in the query dataset

# %%
preds, probs = predict_cell_types_with_rf(adata_query_scvi, adata_ref, "high_level_cell_type")
adata_query_scvi.obs["predicted_high_level_cell_type"] = preds
adata_query_scvi.obs["predicted_high_level_cell_type_probability"] = probs
adata_query_scvi.obs["predicted_high_level_cell_type_probability"].describe()

# %% [markdown]
# Accuracy of the high-level cell type predictions

# %%
# Get the labels for each cell in adata_query_scvi, excluding "unknown" cell types
mask = adata_query_scvi.obs['high_level_cell_type'].astype(str) != "unknown"
original_cell_types = adata_query_scvi.obs['high_level_cell_type'].astype(str)[mask]
predicted_cell_types = adata_query_scvi.obs['predicted_high_level_cell_type'].astype(str)[mask]

# Compute the fraction of correctly predicted cell types (excluding "unknown")
frac_correct = sum(predicted_cell_types == original_cell_types) / len(predicted_cell_types)
print(f"Fraction of correctly predicted high-level cell types (excluding 'unknown'): {frac_correct:.3f}")

# %% [markdown]
# ### Predict low-level cell types for each high-level cell type

# %%
# For each high-level cell type, annotate within that group
col_prev_pred = "predicted_high_level_cell_type"
unique_types = adata_query_scvi.obs[col_prev_pred].dropna().unique()
for ct in unique_types:
    print(f"  Annotating {ct}...")
    # Mask for query and reference cells of this type at previous level
    mask_query = adata_query_scvi.obs[col_prev_pred] == ct
    mask_ref = adata_ref.obs["high_level_cell_type"] == ct
    if mask_ref.sum() == 0 or mask_query.sum() == 0:
        print(f"    Skipping {ct}: no cells in reference or query.")
        continue
    preds, probs = predict_cell_types_with_rf(
        adata_query_scvi, adata_ref, "cell_type", mask_query=mask_query, mask_ref=mask_ref
    )
    # Assign only to the relevant subset
    adata_query_scvi.obs.loc[mask_query, "predicted_cell_type_2"] = preds
    adata_query_scvi.obs.loc[mask_query, "predicted_cell_type_2_probability"] = probs

# %% [markdown]
# Accuracy of the hierarchical predictions

# %%
# Get the labels for each cell in adata_query_scvi, excluding "unknown" cell types
mask = adata_query_scvi.obs['cell_type'].astype(str) != "unknown"
original_cell_types = adata_query_scvi.obs['cell_type'].astype(str)[mask]
predicted_cell_types = adata_query_scvi.obs['predicted_cell_type_2'].astype(str)[mask]

# Compute the fraction of correctly predicted cell types (excluding "unknown")
frac_correct = sum(predicted_cell_types == original_cell_types) / len(predicted_cell_types)
print(f"Fraction of correctly predicted cell types (excluding 'unknown'): {frac_correct:.3f}")

# %% [markdown]
# ## Cluster the test dataset and compute DEGs for input into CASSIA

# %%
adata_query_scvi

# %% [markdown]
# I don't need to worry about QC, as it has already been done by CellxGene: see THYMUBROWN-56. I'll also skip PCA and compute a neighbor graph based on the scvi embedding. This should ensure that the clusters correspond to the previously annotated cell types as closely as possible.

# %%
sc.pp.neighbors(adata_query_scvi, n_neighbors=15, use_rep="scvi2", metric="correlation")

# %%
adata_query_scvi

# %% [markdown]
# Low-resolution clustering to get high-level cell types

# %%
sc.tl.leiden(adata_query_scvi, key_added="leiden_0", resolution=0.1, flavor="igraph", n_iterations=2, random_state=350)
adata_query_scvi.obs[["leiden_0"]].value_counts()

# %% [markdown]
# Compare low-resolution clusters to high-level cell types

# %%
# Assess the overall correspondence between leiden_0 and high_level_cell_type

# Calculate the adjusted Rand index (ARI) and normalized mutual information (NMI)
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# Get the labels for each cell in adata_query_scvi
leiden_labels = adata_query_scvi.obs['leiden_0'].astype(str)
celltype_labels = adata_query_scvi.obs['high_level_cell_type'].astype(str)

# Compute ARI
ari = adjusted_rand_score(celltype_labels, leiden_labels)
print(f"Adjusted Rand Index (ARI) between leiden_0 and high_level_cell_type: {ari:.3f}")

# Compute NMI
nmi = normalized_mutual_info_score(celltype_labels, leiden_labels)
print(f"Normalized Mutual Information (NMI) between leiden_0 and high_level_cell_type: {nmi:.3f}")

# Optionally, display a heatmap of the contingency table for visual assessment
import seaborn as sns
import matplotlib.pyplot as plt

contingency = pd.crosstab(celltype_labels, leiden_labels)
plt.figure(figsize=(10, 6))
sns.heatmap(contingency, annot=True, fmt="d", cmap="viridis")
plt.title("Contingency Table: high_level_cell_type vs leiden_0")
plt.ylabel("high_level_cell_type")
plt.xlabel("leiden_0")
plt.tight_layout()
plt.show()

# %% [markdown]
# The resolution of 0.1 is insufficient, as it fails to separate muscle from connective tissue cells. I'll increase it to 0.2.

# %%
sc.tl.leiden(adata_query_scvi, key_added="leiden_0", resolution=0.2, flavor="igraph", n_iterations=2, random_state=350)
adata_query_scvi.obs[["leiden_0"]].value_counts()

# %% [markdown]
# Compare low-resolution clusters to high-level cell types

# %%
# Assess the overall correspondence between leiden_0 and high_level_cell_type

# Calculate the adjusted Rand index (ARI) and normalized mutual information (NMI)
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# Get the labels for each cell in adata_query_scvi
leiden_labels = adata_query_scvi.obs['leiden_0'].astype(str)
celltype_labels = adata_query_scvi.obs['high_level_cell_type'].astype(str)

# Compute ARI
ari = adjusted_rand_score(celltype_labels, leiden_labels)
print(f"Adjusted Rand Index (ARI) between leiden_0 and high_level_cell_type: {ari:.3f}")

# Compute NMI
nmi = normalized_mutual_info_score(celltype_labels, leiden_labels)
print(f"Normalized Mutual Information (NMI) between leiden_0 and high_level_cell_type: {nmi:.3f}")

# Optionally, display a heatmap of the contingency table for visual assessment
import seaborn as sns
import matplotlib.pyplot as plt

contingency = pd.crosstab(celltype_labels, leiden_labels)
plt.figure(figsize=(10, 6))
sns.heatmap(contingency, annot=True, fmt="d", cmap="viridis")
plt.title("Contingency Table: high_level_cell_type vs leiden_0")
plt.ylabel("high_level_cell_type")
plt.xlabel("leiden_0")
plt.tight_layout()
plt.show()

# %% [markdown]
# At this higher resolution, epithelial cells split into more clusters, but we still don't separate muscle from connective tissue cells. Perhaps that's OK because they are both mesenchymal in origin. Also, some epithelial cells cluster with endothelial, as with the previous resolution.

# %% [markdown]
# We have seen before (THYMUBROWN-56) that Leiden clusters don't exactly correspond to annotated cell types and I can't fully reproduce the Yayon et al. cell type annotation process. That said, the correspondence between clusters and Yayon cell types is pretty good overall. Furthermore, in a new dataset, I wouldn't have access to pre-annotated cell types in the first place. 
# 
# For now, I'll proceed with the 12 high-level clusters obtained at the resolution of 0.2. I'll attempt to use CASSIA's subclustering agent to get low-level cell types within each high-level cluster. Another optional agent to consider would be RAG, which can be used to infer "very detailed or novel annotations". See [CASSIA documentation](https://docs.cassia.bio/en/docs/python/introduction-to-optional-agents) for further details.

# %% [markdown]
# ## Use CASSIA to annotate the clusters

# %% [markdown]
# ### Get a data frame of DE genes

# %% [markdown]
# Normalize the data prior to DEG computation

# %%
# Saving count data
adata_query_scvi.layers["counts"] = adata_query_scvi.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata_query_scvi)
# Logarithmize the data
sc.pp.log1p(adata_query_scvi)
adata_query_scvi

# %% [markdown]
# Compute and plot cluster-specific differentially expressed genes

# %%
sc.tl.rank_genes_groups(adata_query_scvi, groupby="leiden_0", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_query_scvi, groupby="leiden_0", standard_scale="var", n_genes=5, gene_symbols="feature_name")

# %% [markdown]
# Get a data frame of markers for CASSIA

# %%
markers = sc.get.rank_genes_groups_df(adata_query_scvi, group=None)  # Get all groups
print(markers.head())

# %%
markers.describe()

# %% [markdown]
# ### Run CASSIA

# %% [markdown]
# Set and validate API key

# %%
# recommended: read key from environment, do not hardcode
key = os.environ.get("CASSIA_OPENAI_API_KEY")
if not key:
    raise EnvironmentError("Environment variable CASSIA_OPENAI_API_KEY is not set. Please set it before running this cell.")
CASSIA.set_api_key(key, provider="openai")
# then validate
CASSIA.validate_api_keys("openai", force_revalidate=True)

# %% [markdown]
# Run CASSIA in batch mode

# %%
CASSIA.runCASSIA_batch(
    marker=markers,
    output_name="cassia-12-clusters",
    model="best",
    provider="openai",
    reasoning="medium",  # Recommended for balanced speed/quality
    tissue="thymus",
    species="human",
    max_workers=72 # Recommend ~75% of CPU cores
)

# %% [markdown]
# Score the predictions

# %%
CASSIA.runCASSIA_score_batch(
    input_file = "cassia-12-clusters_summary.csv",  # JSON auto-detected
    model = "best",
    provider = "openai",
    max_workers=72
)
results = pd.read_csv("cassia-12-clusters_summary_scored.csv")
results.head()

# %% [markdown]
# For each cluster, print out the following data:
# 1) From adata_query_scvi, a table of high_level_cell_type and cell_type with counts.
# 2) From results, Predicted General Cell Type, Predicted Detailed Cell Type, Possible Mixed Cell Types, and score.

# %%
# For each leiden_0 cluster, print counts of (high_level_cell_type, cell_type)
# and the corresponding CASSIA results (general/detailed/mixed/score) if found.

# Helper to find a results column by keywords (case-insensitive, substring match)
def find_col_by_keywords(df, keywords):
    for col in df.columns:
        lc = col.lower()
        if any(k in lc for k in keywords):
            return col
    return None

# Identify candidate columns in results for the requested fields
col_general = find_col_by_keywords(results, ["predicted general", "general cell", "general_cell", "predicted_general"])
col_detailed = find_col_by_keywords(results, ["predicted detailed", "detailed cell", "detailed_cell", "predicted_detailed"])
col_mixed = find_col_by_keywords(results, ["possible mixed", "mixed cell", "possible_mixed", "mixed"])
col_score = find_col_by_keywords(results, ["score", "cassia_score", "confidence"])

requested_cols = [c for c in (col_general, col_detailed, col_mixed, col_score) if c is not None]

# Try to find which column in results links to the cluster id (common names include 'cluster', 'cluster_id', 'cluster_name', 'leiden')
cluster_link_col = find_col_by_keywords(results, ["cluster", "leiden", "cluster_id", "cluster name", "cluster_name", "clusterid"])

clusters = sorted(adata_query_scvi.obs['leiden_0'].astype(str).unique(), key=lambda x: (int(x) if x.isdigit() else x))

for cl in clusters:
    print(f"\n=== Cluster: {cl} ===")

    # Table from adata_query_scvi
    sub = adata_query_scvi[adata_query_scvi.obs['leiden_0'].astype(str) == cl]
    if sub.n_obs == 0:
        print("No cells for this cluster in adata_query_scvi.")
    else:
        table = sub.obs.groupby(['high_level_cell_type', 'cell_type']).size().reset_index(name='count')
        table = table[table['count'] >= 10].reset_index(drop=True)
        if table.empty:
            print("No high_level_cell_type / cell_type info for cells in this cluster.")
        else:
            print("\nCounts (high_level_cell_type, cell_type):")
            print(table.to_string(index=False))

    # Look up results rows for this cluster
    matched = None
    if cluster_link_col is not None:
        # exact match on the linking column (cast to str)
        mask = results[cluster_link_col].astype(str) == cl
        matched = results[mask]
    if matched is None or matched.empty:
        # try to find rows where any column contains the cluster string
        mask_any = results.apply(lambda row: row.astype(str).str.contains(fr'\b{cl}\b', case=False, regex=True).any(), axis=1)
        matched = results[mask_any]

    if matched is None or matched.empty:
        print("\nNo matching CASSIA row found for this cluster in 'results'. Showing available requested fields (if present) for inspection:")
        if requested_cols:
            print(results[requested_cols].head().to_string(index=False))
        else:
            print("Requested CASSIA columns not found in 'results'. Columns available:")
            print(list(results.columns))
    else:
        print("\nCASSIA predictions for this cluster:")
        if requested_cols:
            for idx, row in matched[requested_cols].iterrows():
                print(f"\nRow {idx}:")
                print(row.to_frame(name="value").to_string())
        else:
            # If we couldn't identify the requested cols, print the whole matched row(s)
            print(matched.to_string(index=False))

# %% [markdown]
# ## A summary table of the results
# 
# By and large, CASSIA aggrees with the CellxGene cell types: 
# | Cluster | CellxGene cell types | CASSIA cell types | CASSIA score |
# |---:|---|---|---|
# | 0 | Mast + B + T + thymocytes | Mast cells (mast cell lineage), strongly proliferating (cycling/progenitor-like) | 85 |
# | 1 | T + thymocytes | αβ T lymphocyte (thymic T cell) | 90 |
# | 2 | plasmacytoid dendritic cells + early T lineage precursor | Thymic plasmacytoid dendritic cells (LAMP5+ pDC subset) | 94 |
# | 3 | B cells | B cells | 92 |
# | 4 | Myeloid cells: macrophages, dendritic etc. | Innate myeloid cells | 90 |
# | 5 | endothelial and some mesothelial | endothelial and mesothelial-like | 90 |
# | 6 | smooth muscle + pericyte | smooth muscle–dominant mural cell population; pericyte–VSMC continuum | 94 |
# | 7 | fibroblast | Fibroblast / PDGFRA+ mesenchymal stromal cell | 93 |
# | 8 | myo/neuro/ciliated mTEC | Thymic myoid + neural | 90 |
# | 9 | Schwann | Schwann | 92 |
# | 10 | cTEC+mcTEC+mTEC | mTEC | 94 |
# | 11 | erythrocyte | Erythroid lineage (nucleated red blood cell / erythroid progenitor) | 95 |


