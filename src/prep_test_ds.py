"""
Add noise and modify cell ids in a test dataset
"""

# %% Import libraries
# Standard python libraries
import numpy as np
import os
import pandas as pd
import warnings

# scVerse
import anndata as ad
from anndata import ImplicitModificationWarning

# Local utilities
from utils import load_config, get_data_dir, initialize_output_directory, save_package_versions

# Some settings to avoid errors/warnings
# - suppress "Transforming to str index" warnings when loading CellxGene census 
#   data into anndata
warnings.filterwarnings("ignore", category=ImplicitModificationWarning)

# - enable writing nullable string arrays to h5ad files
pd.set_option("mode.string_storage", "python")
ad.settings.allow_write_nullable_strings = True
