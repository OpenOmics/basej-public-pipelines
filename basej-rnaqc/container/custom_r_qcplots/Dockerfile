# custom_r_qcplots container (RNA QC version)
# Provides: R 4.5.x with extensive package set (ohchibi, ggpubr, dplyr, etc.)
#           plus Python 3 with pyarrow for Parquet output
# Scripts: rna_qc_plot.R, htseq_summarydf_to_matrix.R
# Tag: custom_r_qcplots_0.2.3
#
# Build:
#   docker build -t basejumper_custom_r_qcplots_0.2.3 .
#
# Note: This is a full rebuild Dockerfile. Build time is significant (~30-45 min)
# due to R package compilation.

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=Etc/UTC

# System update and locale setup
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends locales && \
    locale-gen en_US.UTF-8 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install system dependencies
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        dirmngr \
        gnupg \
        wget \
        curl \
        build-essential \
        git \
        pkg-config \
        ca-certificates \
        # R dependencies
        libcurl4-openssl-dev \
        libuv1-dev \
        libxml2-dev \
        libssl-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libpcre2-dev \
        libcairo2-dev \
        libxt-dev \
        libx11-dev \
        libxpm-dev \
        libxmu-dev \
        libblas-dev \
        liblapack-dev \
        libopenblas-dev \
        libgfortran5 \
        # Python packages
        python3 \
        python3-dev \
        python3-pip \
        python3-setuptools \
        python3-wheel \
        python3-numpy \
        python3-scipy \
        python3-pandas \
        python3-matplotlib \
        python3-h5py \
        python3-seaborn \
        python3-sklearn \
        python3-networkx \
        python3-biopython \
        python3-joblib \
        python3-tqdm \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install pyarrow via pip
RUN pip3 install --break-system-packages --no-cache-dir pyarrow

# Add CRAN repository and install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | gpg --dearmor -o /etc/apt/trusted.gpg.d/cran_ubuntu_key.gpg && \
    echo "deb [signed-by=/etc/apt/trusted.gpg.d/cran_ubuntu_key.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
        > /etc/apt/sources.list.d/r-ubuntu.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror
ENV R_CRAN_MIRROR=https://cloud.r-project.org/

# Install core R packages
RUN R -e "install.packages(c( \
        'devtools', \
        'remotes', \
        'BiocManager', \
        'rlang', \
        'ggplot2', \
        'reshape2', \
        'dplyr', \
        'tidyr', \
        'magrittr', \
        'data.table', \
        'RColorBrewer', \
        'cowplot', \
        'gridExtra', \
        'optparse', \
        'cli', \
        'Hmisc', \
        'ape' \
    ), repos = Sys.getenv('R_CRAN_MIRROR'))"

# Install additional R dependencies
RUN R -e "install.packages(c( \
        'vegan', \
        'scales', \
        'egg', \
        'tibble', \
        'patchwork', \
        'jsonlite', \
        'pheatmap', \
        'ggpubr' \
    ), repos = Sys.getenv('R_CRAN_MIRROR'))"

# Install Bioconductor packages
RUN R -e "BiocManager::install(ask = FALSE, update = FALSE)" && \
    R -e "BiocManager::install('ggtree', ask = FALSE, update = FALSE)"

# Install harrietr (depends on ggtree)
RUN R -e "install.packages('harrietr', repos = Sys.getenv('R_CRAN_MIRROR'), dependencies = TRUE)"

# Install ohchibi from GitHub (pinned to specific commit for reproducibility)
RUN R -e "remotes::install_github('isaisg/ohchibi', \
        ref = '9d3dac92081993508facc58877fdf0cfa8bba77c', \
        dependencies = TRUE, \
        upgrade = 'never', \
        force = TRUE)"

# Install plotting and color palette packages
RUN R -e "install.packages(c('pals', 'gt', 'paletteer'), repos = Sys.getenv('R_CRAN_MIRROR'))"

# Install nloptr/lme4/car/rstatix dependency chain for ggpubr
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends cmake && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    R -e "install.packages('nloptr',  repos = Sys.getenv('R_CRAN_MIRROR'))" && \
    R -e "install.packages('lme4',    repos = Sys.getenv('R_CRAN_MIRROR'))" && \
    R -e "install.packages('car',     repos = Sys.getenv('R_CRAN_MIRROR'))" && \
    R -e "install.packages('rstatix', repos = Sys.getenv('R_CRAN_MIRROR'))" && \
    R -e "install.packages('ggpubr',  repos = Sys.getenv('R_CRAN_MIRROR'))"

# Copy QC plotting scripts
COPY scripts/ /usr/local/bin/

WORKDIR /data

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser && chown appuser:appuser /data
USER appuser

CMD ["/bin/bash"]
