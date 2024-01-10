#!/bin/bash
#
# copy the tar file
cp /staging/bmoore22/cellxgene.tar.gz ./
#
# job exit if any command returns with non-zero exit status (aka failure)
set -e
# set environment name
ENVNAME=cellxgene
export ENVDIR=$ENVNAME
# set up the environment
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate
#
# unzip input files
chmod 777 fine_tuned_geneformer.tar.gz
tar -xzvf fine_tuned_geneformer.tar.gz
mkdir test_data
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz -C test_data/
# Command for myprogram, which will use files from the working directory
python geneformer_example_mapping.py
#
# mkdir for output files
mkdir geneformer_output
mv *.png geneformer_output/
mv *.h5ad geneformer_output/
# tar output files
tar -czvf geneformer_output.tar.gz geneformer_output/
chmod 777 geneformer_output.tar.gz
# Before the script exits, make sure to remove the file
rm cellxgene.tar.gz
rm fine_tuned_geneformer.tar.gz
rm pbmc3k_filtered_gene_bc_matrices.tar.gz
#