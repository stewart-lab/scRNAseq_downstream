# update conda
conda update -n base -c conda-forge conda
# create environment
conda create -n cellxgene python=3.11
# activate
conda activate cellxgene
# install
pip install -U cellxgene-census

# install aws
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
./aws/install -i /w5home/bmoore/aws-cli -b /w5home/bmoore/bin
# test install
/w5home/bmoore/bin/aws --version

# get geneformer
git lfs install
# for chtc follow these instructions: https://gist.github.com/pourmand1376/bc48a407f781d6decae316a5cfa7d8ab
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
# activate conda
conda activate cellxgene
# install geneformer
pip install .
# install torch
pip3 install torch torchvision torchaudio torchdata
pip install transformers[torch]
pip3 install leidenalg
# install pandas and matplot lib
conda install pandas matplotlib

# deactivate
conda deactivate
# pack environment to run on chtc
conda install -c conda-forge conda-pack
conda pack -n cellxgene --dest-prefix='$ENVDIR'
chmod 644 cellxgene.tar.gz
ls -sh cellxgene.tar.gz
# tarball is 3.2 G so move to staging
mv cellxgene.tar.gz /staging/bmoore22/
