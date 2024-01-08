# import packages
import warnings

warnings.filterwarnings("ignore")

from transformers import BertForSequenceClassification
from transformers import Trainer
from geneformer import DataCollatorForCellClassification
from geneformer import TranscriptomeTokenizer
from geneformer import EmbExtractor
from cellxgene_census.experimental import get_embedding
from cellxgene_census.experimental.ml.huggingface import GeneformerTokenizer
import datasets
import json
import os
import scanpy as sc
import numpy as np
import cellxgene_census
import tiledbsoma
# Plotting
import matplotlib
import matplotlib.pyplot as plt

# prep test data: #
## Set the index as the ENSEMBL ID and stores it in the obs column "ensembl_id"
## Add read counts to the obs column "n_counts"
## Add an ID column to be used for joining later in the obs column "joinid"

adata = sc.read_10x_mtx("/test_data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))

h5ad_dir = "/test_data/h5ad/"

if not os.path.exists(h5ad_dir):
    os.makedirs(h5ad_dir)

adata.write(h5ad_dir + "pbmcs.h5ad")

# tokenize the test data using Geneformer’s tokenizer
token_dir = "/test_data/tokenized_data/"

if not os.path.exists(token_dir):
    os.makedirs(token_dir)

tokenizer = TranscriptomeTokenizer(custom_attr_name_dict={"joinid": "joinid"})
tokenizer.tokenize_data(
    data_directory=h5ad_dir,
    output_directory=token_dir,
    output_prefix="pbmc",
    file_format="h5ad",
)

# prepare the data from the model
model_dir = "/fine_tuned_geneformer/"
label_mapping_dict_file = os.path.join(model_dir, "label_to_cell_subclass.json")

with open(label_mapping_dict_file, "r") as fp:
    label_mapping_dict = json.load(fp)

# Loading tokenized data
dataset = datasets.load_from_disk(token_dir + "pbmc.dataset")
dataset
# add column
dataset = dataset.add_column("label", [0] * len(dataset))

# load model and run inference
# 1 GPU recommended
# reload pretrained model
model = BertForSequenceClassification.from_pretrained(model_dir)
# create the trainer
trainer = Trainer(model=model, data_collator=DataCollatorForCellClassification())
# use trainer
predictions = trainer.predict(dataset)

# select most likely class
predicted_label_ids = np.argmax(predictions.predictions, axis=1)
predicted_labels = [label_mapping_dict[str(i)] for i in predicted_label_ids]

# Inspect results
# add predicted labels to anndata
adata.obs["predicted_cell_subclass"] = predicted_labels

# Do normalization, variable genes, scaling, pca, nearest enighbors and umap to visualize results
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# add original labels to compare
sc.tl.leiden(adata)
original_cell_types = [
    "CD4-positive, alpha-beta T cell (1)",
    "CD4-positive, alpha-beta T cell (2)",
    "CD14-positive, monocyte",
    "B cell (1)",
    "CD8-positive, alpha-beta T cell",
    "FCGR3A-positive, monocyte",
    "natural killer cell",
    "dendritic cell",
    "megakaryocyte",
    "B cell (2)",
]
adata.rename_categories("leiden", original_cell_types)

# visualize original labels
with plt.rc_context():
    sc.pl.umap(adata, color="leiden", title="Original Annotations")
    plt.savefig(data_dir + "original_annotation.png")

# visualize predicted labels
with plt.rc_context():
    sc.pl.umap(adata, color="predicted_cell_subclass", title="Predicted Geneformer Annotations")
    plt.savefig(data_dir + "predicted_annotation.png")

