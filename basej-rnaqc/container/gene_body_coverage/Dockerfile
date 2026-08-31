# gene_body_coverage container
# Provides: Python 3 with RSeQC for gene body coverage analysis
# Tag: gene_body_coverage_0.2.1
#
# Build:
#   docker build -t basejumper_gene_body_coverage_0.2.1 .
#
# Uses Ubuntu 22.04 (Python 3.10) for compatibility with RSeQC 4.0.0

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-dev \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libssl-dev \
        ca-certificates \
    && pip3 install --no-cache-dir \
        RSeQC==4.0.0 \
        pandas==1.4.3 \
        seaborn==0.12.0 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy custom scripts
COPY scripts/ /usr/local/bin/

WORKDIR /data

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser && chown appuser:appuser /data
USER appuser

CMD ["/bin/bash"]
