# %% Import modules
from typing import Dict
import numpy as np
import pandas as pd
import anndata as ad

# %% Enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True

# %% A function to create an example AnnData object
def make_example_adata(
    n_obs: int = 50,
    n_vars: int = 100,
    n_pcs: int = 5,
    seed: int = 0
) -> ad.AnnData:
    """
    Build a small AnnData object with:
      - X: normally distributed expression matrix (float32)
      - layers['counts']: integer counts (Poisson)
      - obs: batch and QC-like fields
      - var: gene metadata
      - obsm['X_pca'] and varm['PCs']: random embeddings/loadings
      - uns: small metadata dict
    """
    rng = np.random.default_rng(seed)

    X = rng.normal(loc=0.0, scale=1.0, size=(n_obs, n_vars)).astype(np.float32)
    counts = rng.poisson(lam=5.0, size=(n_obs, n_vars)).astype(np.int32)

    obs = pd.DataFrame(
        {
            "batch": rng.choice(["A", "B"], size=n_obs),
            "percent_mito": rng.uniform(0.0, 0.25, size=n_obs),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    var = pd.DataFrame(
        {
            "gene_symbol": [f"gene_{i}" for i in range(n_vars)],
            "feature_type": ["gene"] * n_vars,
        },
        index=[f"ENSG{i:05d}" for i in range(n_vars)],
    )

    obsm: Dict[str, np.ndarray] = {"X_pca": rng.normal(size=(n_obs, n_pcs)).astype(np.float32)}
    varm: Dict[str, np.ndarray] = {"PCs": rng.normal(size=(n_vars, n_pcs)).astype(np.float32)}
    uns = {"experiment": {"name": "example_adata", "seed": int(seed)}}

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=var,
        layers={"counts": counts},
        obsm=obsm,
        varm=varm,
        uns=uns,
    )

    # simple QC summaries stored in obs
    adata.obs["n_counts"] = adata.layers["counts"].sum(axis=1)
    adata.obs["n_genes_by_counts"] = (adata.layers["counts"] > 0).sum(axis=1)

    return adata

# %% Construct example AnnData and expose it as `adata`
adata = make_example_adata()

# %% Save the AnnData to an .h5ad file
adata.write_h5ad("example_adata.h5ad")