
FROM stewartlab/scrnaseq_downstream3:v2

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

# nichenetr, for src/nichenet.R -- everything else it needs (Seurat,
# SeuratObject, tidyverse, viridis, digest, DiagrammeR) is already in
# scRNAseq_new. Not on CRAN/Bioconductor (confirmed by actually attempting
# BiocManager::install() first -- "package 'nichenetr' is not available for
# Bioconductor version '3.18'"); it's GitHub-only, same as CellChat/autozyme
# below.
#
# nichenetr pulls in ggpubr -> rstatix/car/pbkrtest/doBy/shadowtext/ggiraph as
# transitive deps, whose source builds need system graphics libraries this
# base image doesn't have -- discovered by actually attempting the build,
# which failed on gdtools ("cairo-ft.h: No such file"), systemfonts
# ("ft2build.h: No such file"), and nloptr ("CMake not found") before any of
# Beth's own code even loads.
RUN /opt/conda/bin/conda install -p /scRNAseq_new -c conda-forge -y \
    cairo freetype fontconfig harfbuzz fribidi cmake pkg-config

# Deriv (a doBy -> ggpubr -> nichenetr transitive dep) fails to compile from
# CRAN source against this R build even with the libraries above: its C code
# calls R_ClosureFormals/Rf_allocLang, internal R APIs only made public
# starting in R 4.4 (we're pinned to R 4.3.3 here, matching everything else in
# scRNAseq_new). The conda-forge binary is built against R 4.3 and sidesteps
# the incompatible source build entirely.
RUN /opt/conda/bin/conda install -p /scRNAseq_new -c conda-forge -y r-deriv

RUN . /scRNAseq_new/bin/activate && \
    Rscript -e "remotes::install_github('saeyslab/nichenetr', upgrade = 'never')"

# DiagrammeRsvg + rsvg (would let nichenet.R's sig_network_maps() ->
# export_graph() render NicheNet's ligand-signaling DiagrammeR graphs to PDF)
# -- attempted and reverted. DiagrammeRsvg depends on the V8 R package, whose
# conda-forge binary only exists for r44 (scRNAseq_new is pinned to R 4.3.3,
# shared by ~13 other methods), and compiling Google's V8 engine from CRAN
# source is its own notoriously fragile rabbit hole -- not worth the risk to
# this shared environment for one optional plot type. nichenet.R's
# sig_network_maps() already wraps both export_graph() call sites in
# tryCatch, so this just means those specific signaling-graph PDFs are
# skipped (logged, not fatal) while every other nichenet output still works.
#
# install.packages()/BiocManager::install() don't fail the build on a
# package install error -- they just warn and leave the package missing --
# so verify explicitly and fail the build hard if anything didn't actually
# install (this is what caught the missing-toolchain issue above).
RUN . /scRNAseq_new/bin/activate && Rscript -e "\
pkgs <- c('gprofiler2', 'devEMF', 'rrvgo', 'org.Hs.eg.db', 'org.Dr.eg.db', 'org.Ss.eg.db', 'nichenetr'); \
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
if (length(missing) > 0) stop('Missing required packages: ', paste(missing, collapse = ', '))"

# `cellchat` conda environment, for src/cellchat.R and src/cellchat_interactions.R.
# A separate environment from scRNAseq_new (same reasoning as sccomp2's own
# environment below): CellChat's dependency tree is large and installing it
# into scRNAseq_new -- shared by ~13 other methods -- would risk version
# conflicts with whatever those already need pinned. run_downstream_toolkit.sh
# already dispatches cellchat/cellchat_interactions to `conda run -n cellchat`
# on this assumption; this is what actually creates that environment.
#
# r-base=4.4, not 4.3: tried 4.3 first (for parity with scRNAseq_new) but hit
# two classes of failure building CellChat's dependency tree from CRAN source
# against it -- (1) current conda-forge builds of several of these packages
# (r-seurat, r-matrix, r-reticulate, r-circlize) only exist for r44+ now, no
# r43 build; (2) Deriv's C code calls R_ClosureFormals/Rf_allocLang, R APIs
# only public since R 4.4, so it can't compile against 4.3 at all. Since this
# env shares nothing with scRNAseq_new, there's no actual parity requirement
# -- 4.4 avoids both problems.
RUN /opt/conda/bin/conda create -n cellchat -c conda-forge -c bioconda -y r-base=4.4

# Compiler toolchain into *this* env's prefix -- same reasoning as the
# scRNAseq_new toolchain install above: each conda env's R has its own
# isolated compiler search path, so this has to be repeated per-environment,
# not just once for scRNAseq_new. Still needed even though most R packages
# below come from conda binaries now: CellChat and autozyme themselves are
# GitHub-only and compile their own C++/Rcpp code at install time.
RUN /opt/conda/bin/conda install -n cellchat -c conda-forge -y \
    c-compiler cxx-compiler fortran-compiler make cmake

# Python + umap-learn into the same conda env (conda envs can hold both R
# and Python packages) -- CellChat's functional/structural similarity
# analysis (computeNetSimilarityPairwise -> netEmbedding) calls into Python
# via reticulate for UMAP. Having Python live in this same env/prefix means
# cellchat.R/cellchat_interactions.R can point reticulate at `Sys.which("python")`
# instead of a hardcoded personal virtualenv path.
RUN /opt/conda/bin/conda install -n cellchat -c conda-forge -y python umap-learn

# CRAN/Bioconductor packages needed by cellchat.R and cellchat_interactions.R,
# installed as prebuilt conda-forge/bioconda binaries (matching r-base=4.4)
# rather than compiled from CRAN source via install.packages()/BiocManager --
# discovered the hard way that compiling this same dependency tree from
# source hits missing system graphics/build libraries (cairo, freetype,
# fontconfig, libuv, zlib headers) one at a time; installing as binaries
# sidesteps all of that in one shot instead of chasing each library
# individually. `remotes` is needed for the two GitHub-only installs below
# (CellChat and autozyme are not on CRAN/Bioconductor, so no conda binary
# exists for them).
#
# svglite/ggpubr/BiocNeighbors/systemfonts/Deriv added after CellChat's own
# install_github attempt failed needing them transitively (CellChat's
# DESCRIPTION imports svglite/ggpubr/BiocNeighbors directly, and those in turn
# pull in systemfonts/Deriv, which hit the same missing-system-library and
# R-API-version problems solved above for scRNAseq_new -- installing them as
# binaries here sidesteps a repeat of that whole chase in this env too).
RUN /opt/conda/bin/conda install -n cellchat -c conda-forge -c bioconda -y \
    r-biocmanager r-patchwork r-seurat r-matrix r-ggalluvial \
    r-reticulate r-remotes r-circlize r-jsonlite bioconductor-complexheatmap \
    r-svglite r-ggpubr bioconductor-biocneighbors r-systemfonts r-deriv

# NMF specifically NOT from conda-forge: its conda-forge build is stuck at
# 0.21.0 (no newer build exists there), but CellChat's DESCRIPTION requires
# NMF >= 0.23.0 -- installing the stale binary let NMF "install" fine but
# CellChat's own namespace load then failed with "namespace 'NMF' 0.21.0 is
# being loaded, but >= 0.23.0 is required". CRAN has 0.28, well past
# CellChat's floor, so compile that from source instead (toolchain is already
# installed above; NMF's own deps are plain CRAN packages, not the
# graphics/system-library-heavy kind that caused the earlier failures).
#
# BiocManager::install(), not install.packages(): NMF depends on Biobase,
# which is Bioconductor-only -- plain install.packages() only points at CRAN
# and can't resolve it ("dependency 'Biobase' is not available"). BiocManager
# sets up combined CRAN+Bioconductor repos so both resolve in one call.
RUN /opt/conda/bin/conda run -n cellchat \
    Rscript -e "BiocManager::install('NMF', update = FALSE, ask = FALSE)"

# `conda run -n cellchat`, not `. .../bin/activate`: named conda envs (unlike
# the path-based scRNAseq_new above) don't get a bin/activate script -- this
# is also exactly how run_downstream_toolkit.sh itself invokes this env at
# runtime (`conda run -n cellchat ...`), same as the existing sccomp2 env.
RUN /opt/conda/bin/conda run -n cellchat \
    Rscript -e "remotes::install_github('jinworks/CellChat', upgrade = 'never')"
RUN /opt/conda/bin/conda run -n cellchat \
    Rscript -e "remotes::install_github('ElliotXie/autozyme', subdir = 'autozyme_r', upgrade = 'never')"

# presto (fast Wilcoxon rank-sum test), used by cellchat_interactions.R's
# call to identifyOverExpressedGenes() -- discovered by actually running the
# script, which errored asking for it by name. GitHub-only, same as
# CellChat/autozyme above.
RUN /opt/conda/bin/conda run -n cellchat \
    Rscript -e "remotes::install_github('immunogenomics/presto', upgrade = 'never')"

RUN /opt/conda/bin/conda run -n cellchat Rscript -e "\
pkgs <- c('CellChat', 'patchwork', 'Seurat', 'Matrix', 'NMF', 'ggalluvial', 'reticulate', 'circlize', 'jsonlite', 'ComplexHeatmap', 'autozyme', 'presto'); \
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
if (length(missing) > 0) stop('Missing required packages: ', paste(missing, collapse = ', '))"

RUN rm -rf ./src/ ./data/ ./config.json
COPY src/ ./src/
COPY data/ ./data/
COPY environments/ ./environments/
COPY config.json ./

CMD ["/bin/bash"]