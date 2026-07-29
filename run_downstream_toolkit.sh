#!/bin/bash

# --version <tag>: which stewartlab/scrnaseq_downstream3 image tag to pull
# and run (docker mode only; ignored for conda mode). Defaults to v2 --
# bump this when a new image has been built and pushed via build_push.sh.
# Captured before the parsing loop below consumes "$@" via shift, so the
# provenance record further down can log the invocation exactly as typed.
ORIGINAL_INVOCATION="$0 $*"
IMAGE_VERSION="v2"
IMAGE_VERSION_EXPLICIT=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --version)
      IMAGE_VERSION="$2"
      IMAGE_VERSION_EXPLICIT=true
      shift 2
      ;;
    --version=*)
      IMAGE_VERSION="${1#*=}"
      IMAGE_VERSION_EXPLICIT=true
      shift
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

echo "Step 1: Importing DATA_DIR from config.json"
CONFIG_FILE="./config.json"

METHOD=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['METHOD'])")
echo "METHOD imported as $METHOD"

# cassia needs matplotlib/seaborn baked into cassia_env, which only landed in
# image v1.2.3+ -- default to that tag for this method specifically, unless
# the caller explicitly passed --version (e.g. to test a newer tag).
if [ "$METHOD" == "cassia" ] && [ "$IMAGE_VERSION_EXPLICIT" == "false" ]; then
    IMAGE_VERSION="v1.2.3"
elif [ "$METHOD" != "cassia" ] && [ "$IMAGE_VERSION_EXPLICIT" == "false" ]; then
    IMAGE_VERSION="v2"
fi

DATA_DIR=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['$METHOD']['DATA_DIR'])")
echo "DATA_DIR imported as $DATA_DIR"

# Cellchat parallelism: resolve task_id and file count once, used by both docker and conda branches
if [ "$METHOD" == "cellchat" ]; then
    TASK_ID=$(python -c "import json; v=json.load(open('$CONFIG_FILE'))['cellchat'].get('task_id'); print('' if v is None else int(v))")
    FILELIST_PATH=$(python -c "import json; c=json.load(open('$CONFIG_FILE'))['cellchat']; print(c['DATA_DIR'] + c['FILELIST'])")
    N_FILES=$(awk 'END{print NR}' "$FILELIST_PATH")
fi

echo "Step 1.5: Setting up SHARED_VOLUME and recording run provenance"
SHARED_VOLUME="./shared_volume"
mkdir -p "$SHARED_VOLUME"
chmod 777 "$SHARED_VOLUME"

# Everything needed to answer "what exactly ran here" a year from now:
# script version (git commit -- src/, data/, and config.json are all
# bind-mounted from the working tree, not baked into the Docker image, so
# the image tag alone doesn't pin which analysis code ran), parameters
# (a copy of config.json), and -- for docker runs -- which exact image
# (tag *and* immutable digest, since a tag can be overwritten by a later
# build/push). Deliberately NOT included: per-R-package version pins
# (sessionInfo() output is a separate, script-level concern -- ask for it
# if you want it added) or input-data checksums (expensive for sequencing
# data; DATA_DIR's resolved path below is a cheap-enough proxy).
# One timestamp, computed once, shared by the provenance filename below and
# (via the RUN_TIMESTAMP export further down) each method script's own
# output directory -- so e.g. run_provenance_<ts>.txt and
# output_gprofiler_<ts>/ carry the exact same timestamp and can be matched
# up by filename alone, instead of each independently calling
# Sys.time()/datetime.now() and landing on different digits.
RUN_TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
PROVENANCE_FILE="$SHARED_VOLUME/run_provenance_$RUN_TIMESTAMP.txt"
GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo "unknown (not a git repo?)")
GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
GIT_DIRTY_FILES=$(git status --porcelain 2>/dev/null)
if [ -n "$GIT_DIRTY_FILES" ]; then
  GIT_STATUS="DIRTY -- uncommitted changes present below; this run cannot be exactly reproduced from git history alone
$GIT_DIRTY_FILES"
else
  GIT_STATUS="clean"
fi

{
  echo "=== run_downstream_toolkit.sh provenance ==="
  echo "timestamp: $(date -Iseconds)"
  echo "invocation: $ORIGINAL_INVOCATION"
  echo "METHOD: $METHOD"
  echo "DATA_DIR (resolved): $(realpath "$DATA_DIR" 2>/dev/null || echo "$DATA_DIR")"
  echo ""
  echo "--- scRNAseq_downstream repo state ---"
  echo "git commit: $GIT_COMMIT"
  echo "git branch: $GIT_BRANCH"
  echo "git status: $GIT_STATUS"
  echo ""
  echo "--- config.json (parameters used for this run) ---"
  cat "$CONFIG_FILE" 2>/dev/null || echo "(could not read $CONFIG_FILE)"
} > "$PROVENANCE_FILE"

echo "Provenance recorded to $PROVENANCE_FILE"

echo "Step 2: Docker or conda environment?"
read -p "Do you want to use the docker container or have you installed the conda environment on your computer? Reply y for docker, N for conda [y/N]: " confirm

if [[ "$confirm" =~ ^[Yy]$ ]]; then
  # This script only ever pulls and runs a pre-built image -- it never
  # builds one. Building (compiling rrvgo's dependencies, etc.) takes 20+
  # minutes and only needs to happen once, by a maintainer, when the
  # Dockerfile changes; see build_push.sh for that step.
  echo "Step 2.1: Pulling stewartlab/scrnaseq_downstream3:$IMAGE_VERSION"
  docker pull "stewartlab/scrnaseq_downstream3:$IMAGE_VERSION"

  IMAGE_DIGEST=$(docker inspect --format='{{index .RepoDigests 0}}' "stewartlab/scrnaseq_downstream3:$IMAGE_VERSION" 2>/dev/null || echo "unavailable (image not pulled from a registry?)")
  {
    echo ""
    echo "--- docker image ---"
    echo "run mode: docker"
    echo "image tag: stewartlab/scrnaseq_downstream3:$IMAGE_VERSION"
    echo "image digest: $IMAGE_DIGEST"
  } >> "$PROVENANCE_FILE"

  echo "Step 3: Running Docker container for downstream processing scripts"
  # PROVENANCE_FILE is a host path (./shared_volume/...); translate it to
  # where /shared_volume actually lands inside the container so the R/python
  # scripts below can append their own package-version info to the same
  # provenance record, not a path that doesn't exist in there.
  # -i (not -it): this container runs one dispatched, non-interactive
  # command and exits -- no interactive shell needs a real terminal here.
  # -t requires a TTY, which fails outright in non-interactive contexts
  # (cron, CI, this being run from a script); dropping it only affects
  # whether R's own console/log text gets ANSI color codes, not the
  # analysis output itself (e.g. ggplot's pdf() device is unaffected either
  # way -- it writes color the same regardless of TTY state).
  # Forwarded from the calling shell's own environment, not read from a
  # --env-file here -- so it only reaches the container if you've already
  # exported OPENAI_API_KEY yourself (e.g. in .bashrc, or `export
  # OPENAI_API_KEY=...` before running this script). Only cassia currently
  # falls back to this (when config.json's cassia.openAI_key is blank);
  # harmless no-op for every other method.
  docker run --userns=host -i --rm \
    -e "PROVENANCE_FILE=/shared_volume/$(basename "$PROVENANCE_FILE")" \
    -e "RUN_TIMESTAMP=$RUN_TIMESTAMP" \
    -e "OPENAI_API_KEY=${OPENAI_API_KEY:-}" \
    -v "$(realpath "$DATA_DIR"):/data/input_data:ro" \
    -v "$(realpath "$SHARED_VOLUME"):/shared_volume" \
    -v "$(realpath "$CONFIG_FILE"):/config.json" \
    -v "$(realpath "./src"):/src" \
    -v "$(realpath "./data"):/data" \
    "stewartlab/scrnaseq_downstream3:$IMAGE_VERSION" /bin/bash -c "
        if [ \"$METHOD\" == \"seurat_mapping\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate
            Rscript /src/seurat_mapping.R'
        elif [ \"$METHOD\" == \"seurat_integration\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate
            Rscript /src/seurat_integrate_v5.R'
        elif [ \"$METHOD\" == \"sccomp\" ]; then
            conda run -n sccomp2 /bin/bash -c 'Rscript src/sccomp.R'
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
        elif [ \"$METHOD\" == \"de\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate
            Rscript src/get_DE_genes.R'
        elif [ \"$METHOD\" == \"de_cond\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate
            Rscript src/get_DE_genes_across_cond.R'
        elif [ \"$METHOD\" == \"gprofiler\" ]; then
            /bin/bash -c '. scRNAseq_new/bin/activate
            Rscript src/gprofiler.r'
        elif [ \"$METHOD\" == \"cassia\" ]; then
            conda run -n cassia_env Rscript /src/cassia.R
        elif [ \"$METHOD\" == \"cellchat\" ]; then
            if [ -z \"$TASK_ID\" ]; then
                for i in \$(seq 1 $N_FILES); do
                    conda run -n cellchat Rscript /src/cellchat.R --task_id \$i > /shared_volume/nohup_cellchat_task\${i}.out 2>&1 &
                done
                wait
                echo \"All $N_FILES cellchat jobs complete\"
            else
                conda run -n cellchat Rscript /src/cellchat.R --task_id $TASK_ID
            fi
        else
            echo \"Unknown METHOD: $METHOD\"
            exit 1
        fi
    "

else
    echo "Step 2.1: run mode is conda (local environment, not containerized)"
    {
      echo ""
      echo "--- run environment ---"
      echo "run mode: conda (local environment, not containerized -- no image tag/digest to record)"
    } >> "$PROVENANCE_FILE"
    # No path translation needed here (unlike the docker branch above) --
    # conda mode runs directly on the host, so the R/python scripts below see
    # the same filesystem this script does.
    export PROVENANCE_FILE="$(realpath "$PROVENANCE_FILE")"
    export RUN_TIMESTAMP
    echo "Step 3: Running the script with conda environments"
    if [ "$METHOD" == "seurat_mapping" ]; then
        conda activate scRNAseq_new
        Rscript src/seurat_mapping.R
    elif [ "$METHOD" == "seurat_integration" ]; then
        conda activate scRNAseq_new
        Rscript src/seurat_integrate_v5.R
    elif [ "$METHOD" == "sccomp" ]; then
        conda activate sccomp
        Rscript src/sccomp.R
    elif [ "$METHOD" == "pseudotime" ]; then
        source .venv/bin/activate
        python src/pseudotime.py
    elif [ "$METHOD" == "realtime" ]; then
        conda activate realtime
        python src/realtime.py
    elif [ "$METHOD" == "celltypeGPT" ]; then
        conda activate scRNAseq_new
        Rscript src/CellTypeGPT.R
    elif [ "$METHOD" == "clustifyr" ]; then
        conda activate scRNAseq_new
        Rscript src/clustifyr.R
    elif [ "$METHOD" == "recluster" ]; then
        conda activate scRNAseq_new
        Rscript src/recluster-and-annotate.R
    elif [ "$METHOD" == "featureplots" ]; then
        conda activate scRNAseq_new
        Rscript src/featureplots.R
    elif [ "$METHOD" == "seurat2ann" ]; then
        conda activate scRNAseq_new
        Rscript src/convert_seurat2anndata.R
    elif [ "$METHOD" == "subset_seurat" ]; then
        conda activate scRNAseq_new
        Rscript src/subset_seurat.R
    elif [ "$METHOD" == "phate" ]; then
        conda activate phate_env
        Rscript src/phate.R
    elif [ "$METHOD" == "sctype" ]; then
        conda activate scRNAseq_new
        Rscript src/scType.R
    elif [ "$METHOD" == "de" ]; then
        conda activate scRNAseq_new
        Rscript src/get_DE_genes.R
    elif [ "$METHOD" == "de_cond" ]; then
        conda activate scRNAseq_new
        Rscript src/get_DE_genes_across_cond.R
    elif [ "$METHOD" == "gprofiler" ]; then
        conda activate scRNAseq_new
        Rscript src/gprofiler.r
    elif [ "$METHOD" == "cassia" ]; then
        conda activate cassia_env
        Rscript src/cassia.R
    elif [ "$METHOD" == "cellchat" ]; then
        conda activate cellchat
        mkdir -p logs
        if [ -z "$TASK_ID" ]; then
            echo "task_id is null: launching $N_FILES parallel cellchat jobs"
            for i in $(seq 1 $N_FILES); do
                nohup Rscript src/cellchat.R --task_id $i > logs/nohup_cellchat_task${i}.out 2>&1 &
            done
            echo "All $N_FILES cellchat jobs launched. Monitor with: tail -f logs/nohup_cellchat_task*.out"
        else
            echo "task_id=$TASK_ID: running single cellchat job"
            nohup Rscript src/cellchat.R --task_id $TASK_ID > logs/nohup_cellchat_task${TASK_ID}.out 2>&1 &
        fi
    elif [ "$METHOD" == "get_sample_ds_from_cellxgene" ]; then
        source activate cellxgene_scvi
        python src/get_sample_ds_from_cellxgene.py
    elif [ "$METHOD" == "assign_high_level_cell_types" ]; then
        source activate cellxgene_scvi
        python src/assign_high_level_cell_types.py
    elif [ "$METHOD" == "prep_test_ds" ]; then
        source activate cellxgene_scvi
        python src/prep_test_ds.py
    elif [ "$METHOD" == "cellxgene_scvi" ]; then
        source activate cellxgene_scvi
        python src/cellxgene_scvi.py
    elif [ "$METHOD" == "cluster_adata" ]; then
        source activate cellxgene_scvi
        python src/cluster_adata.py
    else
        echo "Unknown method: $METHOD"
        exit 1
    fi
fi
