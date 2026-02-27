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
