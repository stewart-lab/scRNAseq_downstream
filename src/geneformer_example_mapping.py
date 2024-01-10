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

adata = sc.read_10x_mtx("./test_data/filtered_gene_bc_matrices/hg19/", var_names="gene_ids")
adata.var["ensembl_id"] = adata.var.index
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["joinid"] = list(range(adata.n_obs))

h5ad_dir = "./test_data/h5ad/"

if not os.path.exists(h5ad_dir):
    os.makedirs(h5ad_dir)

adata.write(h5ad_dir + "pbmcs.h5ad")

# tokenize the test data using Geneformer’s tokenizer
token_dir = "./test_data/tokenized_data/"

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
model_dir = "./fine_tuned_geneformer/"
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
    plt.savefig("original_annotation.png")

# visualize predicted labels
with plt.rc_context():
    sc.pl.umap(adata, color="predicted_cell_subclass", title="Predicted Geneformer Annotations")
    plt.savefig("predicted_annotation.png")

# save adata
adata.write("pbmc.h5ad")

# Using the Geneformer fine-tuned model for data projection
# Generating Geneformer embeddings for 10X PBMC 3K data
# get the number of categories (cell subclasses) present in the model

n_classes = len(label_mapping_dict)

# run the EmbExtractor

output_dir = "./test_data/geneformer_embeddings"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

embex = EmbExtractor(
    model_type="CellClassifier",
    num_classes=n_classes,
    max_ncells=None,
    emb_label=["joinid"],
    emb_layer=0,
    forward_batch_size=30,
    nproc=8,
)

embs = embex.extract_embs(
    model_directory=model_dir,
    input_data_file=token_dir + "pbmc.dataset",
    output_directory=output_dir,
    output_prefix="emb",
)

# e-order the embeddings based on "joinid"

embs = embs.sort_values("joinid")
adata.obsm["geneformer"] = embs.drop(columns="joinid").to_numpy()

# visualize embeddings

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep="geneformer")
sc.tl.umap(adata)

with plt.rc_context():
    sc.pl.umap(adata, color="predicted_cell_subclass", title="10X PBMC 3K in Geneformer")
    plt.savefig("geneformer_embedding.png")

# Join Geneformer embeddings from 10X PBMC 3K data with other census data
# get census data
census = cellxgene_census.open_soma(census_version="2023-12-15")

# Some PBMC data from these collections
# 1. https://cellxgene.cziscience.com/collections/c697eaaf-a3be-4251-b036-5f9052179e70
# 2. https://cellxgene.cziscience.com/collections/f2a488bf-782f-4c20-a8e5-cb34d48c1f7e

dataset_ids = ["fa8605cf-f27e-44af-ac2a-476bee4410d3", "3c75a463-6a87-4132-83a8-c3002624394d"]

adata_census = cellxgene_census.get_anndata(
    census=census,
    measurement_name="RNA",
    organism="Homo sapiens",
    obs_value_filter=f"dataset_id in {dataset_ids}",
    obsm_layers=["geneformer"],
)

# select shared genes
adata_census.var_names = adata_census.var["feature_id"]
shared_genes = list(set(adata.var_names) & set(adata_census.var_names))
adata_census = adata_census[:, shared_genes]

# subset to match size of test data
index_subset = np.random.choice(adata_census.n_obs, size=3000, replace=False)
adata_census = adata_census[index_subset, :]

# join census data to test data
adata_census.obs["dataset"] = "Census - " + adata_census.obs["dataset_id"].astype(str)
adata.obs["dataset"] = "10X PBMC 3K"
adata.obs["cell_type"] = "Predicted - " + adata.obs["predicted_cell_subclass"].astype(str)

adata_joined = sc.concat([adata, adata_census], join="outer", label="batch")

# visualize
sc.pp.neighbors(adata_joined, n_neighbors=10, n_pcs=40, use_rep="geneformer")
sc.tl.umap(adata_joined)

with plt.rc_context():
    sc.pl.umap(adata_joined, color="dataset")
    plt.savefig("joined_datasets.png")

with plt.rc_context():
    sc.pl.umap(adata_joined, color="cell_type")
    plt.savefig("joined_datasets_cell_type.png")

# save adata
adata_joined.write("joined_datasets.h5ad")