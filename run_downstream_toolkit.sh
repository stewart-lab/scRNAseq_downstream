#!/bin/bash

echo "Step 1: Importing DATA_DIR from config.json"
CONFIG_FILE="./config.json"

METHOD=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['METHOD'])")
echo "METHOD imported as $METHOD"

DATA_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['$METHOD']['DATA_DIR'])")
echo "DATA_DIR imported as $DATA_DIR"

pwd

echo "Step 2: Running the script"
if [ "$METHOD" == "seurat_mapping" ]; then
    source activate scRNAseq_new
    Rscript src/seurat_mapping.R
elif [ "$METHOD" == "seurat_integration" ]; then
    source activate scRNAseq_new
    Rscript -e "rmarkdown::render('src/seurat_integrate_v5.Rmd')"
elif [ "$METHOD" == "sccomp" ]; then
    source activate sccomp
    Rscript -e "rmarkdown::render('src/ssccomp.Rmd')"
elif [ "$METHOD" == "pseudotime" ]; then
    source .venv/bin/activate
    python src/pseudotime.py
elif [ "$METHOD" == "realtime" ]; then
    source activate realtime
    python src/realtime.py
elif [ "$METHOD" == "celltypeGPT" ]; then
    source activate scRNAseq_best
    Rscript -e "rmarkdown::render('src/CelltypeGPT.Rmd')"
elif [ "$METHOD" == "clustifyr" ]; then
    source activate scRNAseq_best
    Rscript -e "rmarkdown::render('src/clustifyr.Rmd')"
elif [ "$METHOD" == "recluster" ]; then
    source activate scRNAseq_new
    Rscript src/recluster-and-annotate.R
elif [ "$METHOD" == "featureplots" ]; then
    source activate scRNAseq_best
    Rscript src/featureplots.R
elif [ "$METHOD" == "seurat2ann" ]; then
    source activate scRNAseq_best
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
