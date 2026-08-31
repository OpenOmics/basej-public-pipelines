# custom_r_nf_rnaseq container (consolidated 0.9 + 0.10)
# Provides: R 4.x with packages for HTSeq gene annotation parsing,
#           Qualimap QC report, and PCA/heatmap visualization
# Scripts: htseq_to_gene_df.R, parse_qualimap.R, pca_heatmap_qc_salmon_htseq_nometadata.R
# Tag: custom_r_nf_rnaseq_0.10
#
# Build:
#   docker build -t basejumper_custom_r_nf_rnaseq_0.10 .

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=Etc/UTC

# System dependencies
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        locales \
        software-properties-common \
        dirmngr \
        gnupg \
        wget \
        ca-certificates \
        build-essential \
        pkg-config \
        libcurl4-openssl-dev \
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
    && locale-gen en_US.UTF-8 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | gpg --dearmor -o /etc/apt/trusted.gpg.d/cran_ubuntu_key.gpg && \
    echo "deb [signed-by=/etc/apt/trusted.gpg.d/cran_ubuntu_key.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
        > /etc/apt/sources.list.d/r-ubuntu.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENV R_CRAN_MIRROR=https://cloud.r-project.org/

# Install R packages needed for:
# - HTSeq gene annotation (htseq_to_gene_df.R): dplyr, tidyr, data.table, stringr, readr
# - Qualimap parsing (parse_qualimap.R): dplyr, stringr, readr
# - PCA/heatmap (pca_heatmap_qc_salmon_htseq_nometadata.R): ggplot2, pheatmap, reshape2, RColorBrewer, scales
RUN R -e "install.packages(c( \
        'dplyr', \
        'tidyr', \
        'data.table', \
        'stringr', \
        'readr', \
        'optparse', \
        'ggplot2', \
        'reshape2', \
        'pheatmap', \
        'RColorBrewer', \
        'scales', \
        'gridExtra' \
    ), repos = Sys.getenv('R_CRAN_MIRROR'))"

# Copy scripts
COPY scripts/ /usr/local/bin/

WORKDIR /data

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser && chown appuser:appuser /data
USER appuser

CMD ["/bin/bash"]
