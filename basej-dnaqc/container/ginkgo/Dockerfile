# ginkgo container (Ginkgo CNV: binUnsorted, CNVcaller, R scripts)
# Tag: ginkgo_0.3.1
# This single image bundles every Ginkgo script (binUnsorted, CNVcaller,
# cnv_ginkgo.R, rds_to_flat.R, rscript_parse_rds_cnv_metrics.R) and serves all
# five Ginkgo processes (GINKGO_BINUNSORT, GINKGO_SEGMENTATION_R,
# GINKGO_CNV_CALLER, GINKO_RDS_TO_FLAT, PARSE_RDS_CNV_METRICS).
# Build:
#   docker build -t basejumper_ginkgo_0.3.1 .
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive
COPY ./binUnsorted /usr/local/bin
COPY ./CNVcaller /usr/local/bin

# The COPYed binaries inherit the source files' permissions, which lack the
# execute bit (-rw-rw-r--). binUnsorted and CNVcaller are invoked directly by
# the Ginkgo processes, so without +x they fail with exit 126 ("Permission
# denied"), silently emptying the Ginkgo output channels. Make them executable.
RUN chmod +x /usr/local/bin/binUnsorted /usr/local/bin/CNVcaller


RUN apt-get update -qq -y \
    && apt-get install -y --no-install-recommends software-properties-common dirmngr gnupg
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN add-apt-repository -y  "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/"
RUN apt-get install -y --no-install-recommends r-base r-base-dev build-essential

# Pin CRAN to a dated snapshot from the R 3.6 era (just before R 4.0.0, 2020-04-24)
# so install.packages() resolves to versions compatible with R 3.6.3 instead of
# today's CRAN (which only serves newer-R builds). Uses Posit Public Package
# Manager binary packages for Ubuntu bionic; the HTTPUserAgent line is required
# for PPM to serve binaries (falls back to source compile otherwise).
RUN mkdir -p /usr/lib/R/etc && \
    printf '%s\n' \
      'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' \
      'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/bionic/2020-04-20"))' \
      >> /usr/lib/R/etc/Rprofile.site

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('ctc',version = '3.10')"
RUN R -e "BiocManager::install('DNAcopy',version = '3.10')"
RUN R -e "install.packages('inline')"
RUN R -e "install.packages('fastcluster')"
RUN R -e "install.packages('rlang')"
RUN R -e "install.packages('Rcpp')"
RUN R -e "install.packages('colorspace')"
RUN R -e "install.packages('scales')"
RUN R -e "install.packages('glue')"
RUN R -e "install.packages('ellipsis')"
RUN R -e "install.packages('magrittr')"
RUN R -e "install.packages('fansi')"
RUN R -e "install.packages('utf8')"
RUN R -e "install.packages('vctrs')"
RUN R -e "install.packages('tibble')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('gridExtra')"
RUN R -e "install.packages('heatmap3')"
RUN R -e "install.packages('gtools')"
RUN R -e "install.packages('gplots')"
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('farver')"
RUN R -e "install.packages('digest')"
RUN R -e "install.packages('dplyr')"

COPY ./cnv_ginkgo.R /usr/local/bin
COPY ./rds_to_flat.R /usr/local/bin
COPY ./rscript_parse_rds_cnv_metrics.R /usr/local/bin

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser
USER appuser

CMD ["/bin/bash"]
