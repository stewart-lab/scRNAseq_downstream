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
  echo "Step 2.1: Removing and recreating the SHARED_VOLUME"
  TIMESTAMP=$(date +%Y%m%d_%H%M%S)
  SHARED_VOLUME="./shared_volume_$TIMESTAMP"

  rm -rf "$SHARED_VOLUME"
  mkdir -p "$SHARED_VOLUME"
  chmod 777 "$SHARED_VOLUME"

  echo "Step 2.2: Building Docker container for downstream processing"
  # Build the Docker image from the pre_pipeline directory
  docker build -t scaligner_v2_with_genomes_and_jq ./pre_pipeline

echo "Step 3: Running the script"
if [ "$METHOD" == "seurat_mapping" ]; then
    source activate scRNAseq_new2
    Rscript src/seurat_mapping.R
elif [ "$METHOD" == "seurat_integration" ]; then
    source activate scRNAseq_new2
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
    source activate scRNAseq_best
    Rscript -e "rmarkdown::render('src/subset_seurat.Rmd')"
elif [ "$METHOD" == "phate" ]; then
    source activate phate_env
    Rscript -e "rmarkdown::render('src/phate.Rmd')"
elif [ "$METHOD" == "sctype" ]; then
    source activate scRNAseq_new
    Rscript src/scType.R
else
    echo "Unknown method: $METHOD"
    exit 1
fi
