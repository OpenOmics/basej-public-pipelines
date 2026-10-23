# custom_snp_somatic_filter_sequoia container (feb2026)
# Provides: R 4.2.0 stack (sequoia/treemut/ggtree phylogeny, data.table, VGAM,
#           GenomicRanges/Rsamtools), bcftools, gawk, poppler-utils, python3 +
#           MuSiCaL, and MutationalPatterns (SigDyn). Bundles every R/awk/py
#           script used by the somatic filtering, lineage phylogeny and
#           mutational-signature processes (installed to /usr/local/bin).
# Tag: custom_snp_somatic_filter_sequoia_feb2026
# Source: nf-bioskryb-utils containers/snp_somatic_filter_sequoia
#         (feature/isg/snp_indel_somatic_pipeline_refactoring_invert_order_tests_both_position_withbulkremoval)
#
# Build (x86/amd64 only):
#   docker build -t basejumper_custom_snp_somatic_filter_sequoia_feb2026 .

FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update -qq -y 
RUN apt install -y  --no-install-recommends software-properties-common dirmngr gnupg
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN add-apt-repository -y  "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/"
RUN apt-get update
# Toolchain + headers for compiling CRAN/Bioc from source (same class of deps as Posit/Ubuntu R docs).
# Includes cmake (fs/libuv), gfortran (many stats/BioC packages), zlib/bz2/lzma/pcre (common native deps).
RUN apt-get install -y --no-install-recommends \
    autoconf automake \
    build-essential \
    cmake \
    gfortran \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libicu-dev \
    libjpeg-dev \
    liblzma-dev \
    libpcre2-dev \
    libpng-dev \
    libssh2-1-dev \
    libssl-dev \
    libtiff5-dev \
    libuv1-dev \
    libxml2-dev \
    pkg-config \
    wget \
    zlib1g-dev

# Do not `rm -rf /var/lib/apt/lists/*` here: the next layers run `apt-get install` and need the index.
ENV R_VERSION 4.2.0
RUN apt-get update -qq && apt-get install -y --no-install-recommends gdebi-core
RUN wget  https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb --no-check-certificate
RUN gdebi --non-interactive r-${R_VERSION}_1_amd64.deb
RUN ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
RUN ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install pinned core CRAN stack with remotes.
RUN R -e "install.packages('remotes', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('BiocManager', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('rlang', version = '1.0.2', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('Rcpp', version = '1.0.8.3', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('scales', version = '1.2.0', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('magrittr', version = '2.0.3', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('tibble', version = '3.1.7', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('ggplot2', version = '3.5.2', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('gridExtra', version = '2.3', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('Rmisc', version = '1.5.1', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('plotrix', version = '3.8-2', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('ggrepel', version = '0.9.1', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('dplyr', version = '1.0.9', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('cli', version = '3.3.0', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('readr', version = '2.1.2', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('reshape2', version = '1.4.4', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('ape', version = '5.8', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('seqinr', version = '4.2-36', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('stringr', version = '1.5.1', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('data.table', version = '1.15.4', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('tidyr', version = '1.3.1', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('VGAM', version = '1.1-11', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('MASS', version = '7.3-9', repos = 'https://cloud.r-project.org')"

RUN R -e  "BiocManager::install('Rsamtools',version = '3.16')" 
RUN R -e  "BiocManager::install('GenomicRanges',version = '3.16')" 
# tidytree 0.4.7 dropped random_ref; treeio 1.22 zzz.R still expects it — pin 0.4.6.
# ggplot2 3.5.2 + aplot need exported is_ggplot; ggtree 3.6.2 still used getFromNamespace("rect_to_poly","ggplot2")
# but ggplot2 >=3.5 removed it — patch ggtree geom_hilight.R with a local rect_to_poly (ggplot2 3.4.x implementation).
# ggplot2 3.5 also removed check_linewidth from the namespace — patch with local .check_linewidth (see patch_ggtree_geom_hilight.R).
# ggplot2 3.5 empty() uses is_waiver internally; ggtree 3.6.2 used getFromNamespace("empty") — patch geom_segment.R (YuLab master).
# Install deps, patch, and R CMD INSTALL in one RUN so ggtree sees the same library as Bioc/CRAN installs.
# ggtree Imports include ggfun + yulab.utils (not pulled by treeio alone). Use wget (already installed) for tarball.
COPY patch_ggtree_geom_hilight.R /tmp/patch_ggtree_geom_hilight.R
COPY patch_ggtree_geom_segment.R /tmp/patch_ggtree_geom_segment.R
RUN R -e "remotes::install_version('tidytree', version = '0.4.6', repos = 'https://cloud.r-project.org', upgrade = 'always')" && \
    R -e "install.packages('fs', repos = 'https://cloud.r-project.org')" && \
    R -e "install.packages(c('aplot', 'ggfun', 'yulab.utils'), repos = 'https://cloud.r-project.org')" && \
    R -e "BiocManager::install('treeio', version = '3.16', update = FALSE, ask = FALSE)" && \
    wget -q -O /tmp/ggtree_3.6.2.tar.gz https://bioconductor.org/packages/3.16/bioc/src/contrib/ggtree_3.6.2.tar.gz && \
    cd /tmp && tar -xzf ggtree_3.6.2.tar.gz && \
    perl -i -pe 's/rect_to_poly <- getFromNamespace\("rect_to_poly", "ggplot2"\)/rect_to_poly <- function(xmin, xmax, ymin, ymax) { data.frame(y = c(ymax, ymax, ymin, ymin, ymax), x = c(xmin, xmax, xmax, xmin, xmin), stringsAsFactors = FALSE) }/' /tmp/ggtree/R/geom_hilight.R && \
    Rscript /tmp/patch_ggtree_geom_hilight.R && \
    Rscript /tmp/patch_ggtree_geom_segment.R && \
    R CMD INSTALL /tmp/ggtree
RUN R -e "library(ggtree); library(treeio); cat('ggtree', as.character(packageVersion('ggtree')), 'treeio', as.character(packageVersion('treeio')), 'tidytree', as.character(packageVersion('tidytree')), 'OK\n')"

RUN R -e  "install.packages('optparse', repos = 'http://cran.us.r-project.org')" 
RUN R -e  "install.packages('getopt', repos = 'http://cran.us.r-project.org')" 

RUN apt-get install -y git
# devtools>=2.4 pulls testthat 3 -> waldo; PPM's diffobj for 2022-06-15 is <0.3.5 so waldo 0.4.0 breaks
# (compare_proxy). Pin testthat 2.3 + diffobj 0.3.5, devtools 2.3.2, and dependencies=NA (no Suggests) so
# lintr/mockery do not reinstall testthat 3.
RUN R -e "ppm <- 'https://packagemanager.posit.co/cran/2022-06-15'; options(repos = c(CRAN = ppm)); remotes::install_version('testthat', version = '2.3.2', repos = ppm, upgrade = 'never'); remotes::install_version('diffobj', version = '0.3.5', repos = ppm, upgrade = 'always'); remotes::install_version('devtools', version = '2.3.2', repos = ppm, upgrade = 'never', dependencies = NA); if (!requireNamespace('devtools', quietly = TRUE)) stop('devtools install failed'); devtools::install_git('https://github.com/NickWilliamsSanger/treemut')"

RUN apt-get install bc

RUN R -e  "install.packages('foreach', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('doParallel', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('parallelly', repos = 'http://cran.us.r-project.org')"

RUN R -e  "install.packages('pals', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('paletteer', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('egg', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('RColorBrewer', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('jsonlite', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('png', repos = 'http://cran.us.r-project.org')"
RUN R -e  "install.packages('hexbin', repos = 'http://cran.us.r-project.org')"

# Later CRAN installs can pull tidytree 0.4.7 and break treeio lazy-load; keep 0.4.6.
# Pin ggfun to current CRAN (ggplot2 3.5–aware helpers); then re-pin ggplot2/tidytree so nothing drifts.
RUN R -e "remotes::install_version('ggfun', version = '0.2.0', repos = 'https://cloud.r-project.org', upgrade = 'always')" && \
    R -e "remotes::install_version('tidytree', version = '0.4.6', repos = 'https://cloud.r-project.org', upgrade = 'always')" && \
    R -e "remotes::install_version('ggplot2', version = '3.5.2', repos = 'https://cloud.r-project.org', upgrade = 'always')" && \
    R -e "library(ggtree); library(treeio); cat('ggtree stack re-check OK\n')"

RUN apt install -y  gawk

# SUBSET_ANNOTATED_VCFS_FOR_MUTSIG / EXPORT_MUTSIG_SOMATIC_VARIANT_TABLE run bcftools + awk from this image.
RUN apt-get update -qq && apt-get install -y --no-install-recommends bcftools python3

RUN apt-get install -y poppler-utils

# ── MuSiCaL (Park lab) ────────────────────────────────────────────────────────
# Python package for COSMIC SBS signature refitting (likelihood-bidirectional).
# pip3 is installed alongside the package; MuSiCaL pulls numpy/scipy/pandas/sklearn.
RUN apt-get update -qq && apt-get install -y --no-install-recommends python3-pip
# pip install MuSiCaL resolves to an unrelated music-theory library (v0.2.0).
# The Park lab's bioinformatics tool (parklab/MuSiCaL) must be installed from GitHub.
RUN pip3 install --no-cache-dir numpy scipy pandas scikit-learn matplotlib seaborn tqdm git+https://github.com/parklab/MuSiCaL.git

# ── SigDyn (Goncalves lab) — uses MutationalPatterns for COSMIC fitting ───────
# MutationalPatterns provides fit_to_signatures_strict() and get_known_signatures()
# with bundled COSMIC v3.3.1 SBS signatures (no internet needed at run time).
RUN R -e "BiocManager::install('MutationalPatterns', version = '3.16', update = FALSE, ask = FALSE)"

# ── musicatk (Campbell lab) — NOT INSTALLED ───────────────────────────────────
# musicatk requires R >= 4.4.0 (GitHub HEAD) and its BioConductor 3.16 release
# has unresolvable deps (topicmodels compile failure, FactoMineR unavailable for
# R 4.2).  This container uses R 4.2.0.  musicatk_enabled must be false when
# using this image; a separate R 4.4+ image is needed to run musicatk.

# Export COSMIC v3.2 SBS signatures from MutationalPatterns (already installed,
# no network needed).  MutationalPatterns 3.8.x (BioConductor 3.16) supports up
# to COSMIC_v3.2; v3.3+ is not available in this Bioc version.
# Consumed by run_musical_analysis.py (Python) and run_sigdyn_analysis.R (R).
RUN mkdir -p /usr/local/share/cosmic && \
    Rscript -e "library(MutationalPatterns); s <- get_known_signatures(source='COSMIC_v3.2'); write.table(s, '/usr/local/share/cosmic/COSMIC_v3.2_SBS_WGS.tsv', sep='\t', quote=FALSE, col.names=NA)" && \
    ls -lh /usr/local/share/cosmic/

# Verify the three new tools load cleanly
RUN python3 -c "\
import musical, sys; \
ver = getattr(musical, '__version__', '?'); \
has_refit = hasattr(musical, 'refit'); \
has_lc = hasattr(musical, 'load_catalog'); \
print(f'MuSiCaL {ver}  refit={has_refit}  load_catalog={has_lc}'); \
assert has_refit, 'musical.refit missing — wrong package installed'; \
print('MuSiCaL OK') \
"
RUN R -e "library(MutationalPatterns); cat('MutationalPatterns', as.character(packageVersion('MutationalPatterns')), 'OK\n')"
RUN ls -lh /usr/local/share/cosmic/

# Entire scripts/ tree is installed to PATH (/usr/local/bin): R scripts, mutsig *.awk, run_matrix_generator.py, etc.
COPY ./scripts /usr/local/bin

RUN chmod +x /usr/local/bin/run_matrix_generator.py || true
RUN chmod +x /usr/local/bin/mpboot-avx


CMD ["/bin/bash"]
