
FROM stewartlab/scrnaseq_downstream3:v1

# Docker's default RUN shell (sh/dash) confuses conda's activate script's
# shell-detection ("Unrecognized shell."); switch to bash for the
# `. /scRNAseq_new/bin/activate` steps below.
SHELL ["/bin/bash", "-c"]

# scRNAseq_new's R is a conda-forge build, so its Makeconf hardcodes
# conda's own cross-compiler names (x86_64-conda-linux-gnu-cc/-c++/
# -gfortran) rather than plain system gcc/g++ -- a system apt
# build-essential install (tried first, then reverted) provides /usr/bin/
# gcc etc., which R's Makeconf never looks for, so it doesn't help. Any R
# package needing native compilation -- several of rrvgo's own
# dependencies (slam, wordcloud, GOSemSim, tm, umap) among them -- fails
# silently either way (install.packages()/BiocManager::install() don't
# fail the build; they just emit a warning and leave the package
# missing) until the matching conda-forge compiler packages are installed
# into the *same* prefix R itself lives in.
RUN /opt/conda/bin/conda install -p /scRNAseq_new -c conda-forge -y \
    c-compiler cxx-compiler fortran-compiler make

# Packages needed by src/gprofiler.r that aren't in the base image's
# scRNAseq_new environment: gprofiler2/devEMF are needed by the script's
# existing GO-enrichment step (it was already broken in Docker without
# these), and rrvgo + OrgDb annotation packages support the GO-term
# redundancy reduction step. Installed into scRNAseq_new specifically,
# since that's the R environment run_downstream_toolkit.sh activates for
# gprofiler.r. OrgDb packages installed here must stay in sync with
# data/organism_orgdb_map.txt -- add both a package here and a mapping row
# there when supporting a new organism.
#
# Must `. /scRNAseq_new/bin/activate` (not just invoke
# /scRNAseq_new/bin/Rscript by absolute path) so PATH includes
# /scRNAseq_new/bin: R itself runs fine either way, but when it shells out
# to run make/the compiler for a source package, that subprocess's PATH
# needs to contain the directory those binaries actually live in, or
# compilation fails with "make: not found" despite make genuinely existing
# on disk.
RUN . /scRNAseq_new/bin/activate && \
    Rscript -e "install.packages(c('gprofiler2', 'devEMF'), repos = 'https://cloud.r-project.org')"
RUN . /scRNAseq_new/bin/activate && \
    Rscript -e "BiocManager::install(c('rrvgo', 'org.Hs.eg.db', 'org.Dr.eg.db', 'org.Ss.eg.db'), update = FALSE, ask = FALSE)"

# install.packages()/BiocManager::install() don't fail the build on a
# package install error -- they just warn and leave the package missing --
# so verify explicitly and fail the build hard if anything didn't actually
# install (this is what caught the missing-toolchain issue above).
RUN . /scRNAseq_new/bin/activate && Rscript -e "\
pkgs <- c('gprofiler2', 'devEMF', 'rrvgo', 'org.Hs.eg.db', 'org.Dr.eg.db', 'org.Ss.eg.db'); \
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
if (length(missing) > 0) stop('Missing required packages: ', paste(missing, collapse = ', '))"

RUN rm -rf ./src/ ./data/ ./config.json
COPY src/ ./src/
COPY data/ ./data/
COPY environments/ ./environments/
COPY config.json ./

CMD ["/bin/bash"]