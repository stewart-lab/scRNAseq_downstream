#!/bin/bash
# Build and publish the scRNAseq_downstream toolkit image.
#
# Run this by hand whenever the Dockerfile (or a script/package it installs)
# changes -- not on every analysis run. run_downstream_toolkit.sh only pulls
# and runs a tagged image; it never builds. Bump the version below when
# publishing a new image so a specific analysis run stays reproducible
# (":latest"/a floating tag would mean the same command could behave
# differently after a later rebuild).
set -euo pipefail

IMAGE="stewartlab/scrnaseq_downstream3"
VERSION="v2"

echo "Building ${IMAGE}:${VERSION} from Dockerfile..."
docker build -t "${IMAGE}:${VERSION}" .

echo "Pushing ${IMAGE}:${VERSION} to Docker Hub..."
docker push "${IMAGE}:${VERSION}"

echo "Done. Update run_downstream_toolkit.sh's image tag to ${VERSION} if this is the new default."
