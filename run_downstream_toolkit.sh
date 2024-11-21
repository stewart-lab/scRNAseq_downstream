#!/bin/bash

echo "Step 1: Importing DATA_DIR from config.json"
CONFIG_FILE="./config.json"

METHOD=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['METHOD'])")
echo "METHOD imported as $METHOD"

DATA_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['$METHOD']['DATA_DIR'])")
echo "DATA_DIR imported as $DATA_DIR"

echo "Step 2: Docker or conda environment?"
read -p "Do you want to use the docker container or have you installed the conda environment on your computer? Reply y for docker, N for conda [y/N]: " confirm

if [[ "$confirm" =~ ^[Yy]$ ]]; then
  echo "Step 2.1: creating the SHARED_VOLUME"
  SHARED_VOLUME="./shared_volume"

  #rm -rf "$SHARED_VOLUME"
  mkdir -p "$SHARED_VOLUME"
  chmod 777 "$SHARED_VOLUME"

  echo "Step 2.2: Building Docker container for downstream processing"
  # Build the Docker image from the pre_pipeline directory
  docker build -t stewartlab/scrnaseq_downstream2:v1 ./

  echo "Step 3: Running Docker container for downstream processing scripts"
  docker run --userns=host -it \
    -v "$(realpath "$DATA_DIR"):/data/input_data:ro" \
    -v "$(realpath "$SHARED_VOLUME"):/shared_volume" \
    -v "$(realpath "$CONFIG_FILE"):/config.json:ro" \
    scrnaseq_downstream2:v1 /bin/bash -c "
        if [ \"$METHOD\" == \"seurat_mapping\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript /src/seurat_mapping.R'
        elif [ \"$METHOD\" == \"seurat_integration\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript /src/seurat_integrate_v5.R'
        elif [ \"$METHOD\" == \"sccomp\" ]; then
            /bin/bash -c '. sccomp/bin/activate 
            Rscript src/sccomp.R'
        elif [ \"$METHOD\" == \"pseudotime\" ]; then
            /bin/bash -c 'source pst_env/bin/activate
            python src/pseudotime.py'
        elif [ \"$METHOD\" == \"realtime\" ]; then
            /bin/bash -c '. realtime/bin/activate 
            python src/realtime.py'
        elif [ \"$METHOD\" == \"celltypeGPT\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript src/CellTypeGPT.R'
        elif [ \"$METHOD\" == \"clustifyr\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript src/clustifyr.R'
        elif [ \"$METHOD\" == \"recluster\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript src/recluster-and-annotate.R'
        elif [ \"$METHOD\" == \"featureplots\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript src/featureplots.R'
        elif [ \"$METHOD\" == \"seurat2ann\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript src/convert_seurat2anndata.R'
        elif [ \"$METHOD\" == \"subset_seurat\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript src/subset_seurat.R'
        elif [ \"$METHOD\" == \"phate\" ]; then
            conda run -n phate /bin/bash -c 'Rscript src/phate.R'
        elif [ \"$METHOD\" == \"sctype\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate 
            Rscript src/scType.R'
        else
            echo \"Unknown METHOD: $METHOD\"
            exit 1
        fi
    "

else
    echo "Step 2.1: Removing and recreating the SHARED_VOLUME"
    SHARED_VOLUME="./shared_volume"
    #rm -rf "$SHARED_VOLUME"
    mkdir -p "$SHARED_VOLUME"
    chmod 777 "$SHARED_VOLUME"
    echo "Step 3: Running the script with conda environments"
    if [ "$METHOD" == "seurat_mapping" ]; then
        source activate scRNAseq_new
        Rscript src/seurat_mapping.R
    elif [ "$METHOD" == "seurat_integration" ]; then
        source activate scRNAseq_new
        Rscript src/seurat_integrate_v5.R
    elif [ "$METHOD" == "sccomp" ]; then
        source activate sccomp
        Rscript src/sccomp.R
    elif [ "$METHOD" == "pseudotime" ]; then
        source .venv/bin/activate
        python src/pseudotime.py
    elif [ "$METHOD" == "realtime" ]; then
        source activate realtime
        python src/realtime.py
    elif [ "$METHOD" == "celltypeGPT" ]; then
        source activate scRNAseq_new
        Rscript src/CellTypeGPT.R
    elif [ "$METHOD" == "clustifyr" ]; then
        source activate scRNAseq_new
        Rscript src/clustifyr.R
    elif [ "$METHOD" == "recluster" ]; then
        source activate scRNAseq_new
        Rscript src/recluster-and-annotate.R
    elif [ "$METHOD" == "featureplots" ]; then
        source activate scRNAseq_new
        Rscript src/featureplots.R
    elif [ "$METHOD" == "seurat2ann" ]; then
        source activate scRNAseq_new
        Rscript src/convert_seurat2anndata.R
    elif [ "$METHOD" == "subset_seurat" ]; then
        source activate scRNAseq_new
        Rscript src/subset_seurat.R
    elif [ "$METHOD" == "phate" ]; then
        source activate phate_env
        Rscript src/phate.R
    elif [ "$METHOD" == "sctype" ]; then
        source activate scRNAseq_new
        Rscript src/scType.R
    else
        echo "Unknown method: $METHOD"
        exit 1
    fi
fi
