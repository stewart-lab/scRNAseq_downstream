""" 
Utility functions
"""
from datetime import datetime
import json
import os
import shutil
import subprocess

# %% 
# ## Utility functions
# %% Load configuration from config.json file
def load_config():
    """
    Load configuration from config.json file.
    
    Returns
    -------
    dict
        Configuration dictionary loaded from config.json
    """
    work_dir = os.getcwd()
    with open(work_dir + '/config.json') as f:
        config_dict = json.load(f)
        method = config_dict["METHOD"]
        print(
            "loaded parameters from config file: ", 
            config_dict[method]
        )
    return config_dict

# %% Set input data directory based on config settings
# %%
def get_data_dir(config_dict):
    """
    Get the data directory path from configuration.
    
    Parameters
    ----------
    config_dict : dict
        Configuration dictionary containing docker and DATA_DIR settings
    
    Returns
    -------
    str
        Path to the data directory
    """
    docker = config_dict["docker"]
    method = config_dict["METHOD"]
    if docker == "TRUE" or docker == "true" or docker == "T" or docker == "t":
        data_dir = "./data/input_data/"
    else:
        data_dir = config_dict[method]["DATA_DIR"]
    
    return data_dir

# %% Initialize the output directory based on the method name and current timestamp.
def initialize_output_directory(config_dict):
    """
    Initialize the output directory based on the method name and current timestamp.
    
    Parameters
    ----------
    config_dict : dict
        Configuration dictionary
    
    Returns
    -------
    str
        Path to the initialized output directory
    """
    now = datetime.now()
    now = now.strftime("%Y%m%d_%H%M%S")
    out_dir = "./shared_volume/" + config_dict["METHOD"] + "_" + now +"/"
    print("out_dir: ", out_dir)
    os.makedirs(out_dir, mode=0o777, exist_ok=True)
    # copy config file
    shutil.copy('./config.json', out_dir) 

    return out_dir

# %% Save package versions to the output directory
def save_package_versions(out_dir):
    """
    Save package versions from conda and pip to the output directory.
    
    Parameters
    ----------
    out_dir : str
        The output directory where version files will be saved.
    """
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

# %% Compare two cell metadata columns
def compare_cell_metadata_cols(metadata_col1, metadata_col2, adata, out_dir):
    """
    Compare two cell metadata columns: compute similarity metrics and plot the contingency table.
    
    Parameters
    ----------
    metadata_col1 : string
        The name of the first cell metadata column to compare.
    metadata_col2 : string
        The name of the second cell metadata column to compare.
    adata : AnnData
        Annotated data matrix.
    out_dir : str
        Output directory to save the contingency table plot.
    
    Returns
    -------
    tuple[float, float]
        A tuple containing:
        - ari: Adjusted Rand Index
        - nmi: Normalized Mutual Information
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    # Extract the two metadata columns
    col1 = adata.obs[metadata_col1]
    col2 = adata.obs[metadata_col2]

    # Compute and print similarity metrics
    ari = adjusted_rand_score(col1, col2)
    nmi = normalized_mutual_info_score(col1, col2)
    print()
    print(f"Comparing '{metadata_col1}' and '{metadata_col2}':")
    print(f"Adjusted Rand Index (ARI): {ari:.4f}")
    print(f"Normalized Mutual Information (NMI): {nmi:.4f}")

    # Create a contingency table
    contingency_table = pd.crosstab(col1, col2)

    # Plot the contingency table
    plt.figure(figsize=(10, 8))
    sns.heatmap(contingency_table, annot=True, fmt='d', cmap='viridis')
    plt.title(f'Contingency Table: {metadata_col1} vs {metadata_col2}\nARI: {ari:.4f}, NMI: {nmi:.4f}')
    plt.xlabel(metadata_col2)
    plt.ylabel(metadata_col1)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(out_dir, f'{metadata_col1}_vs_{metadata_col2}_contingency.png'), bbox_inches='tight')
    plt.close()

    return ari, nmi