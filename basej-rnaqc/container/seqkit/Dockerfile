# SeqKit (linux/amd64 only).
# Tag: seqkit-2.13.0
# Build:
#   docker build --platform linux/amd64 -t basejumper_seqkit-2.13.0 .
FROM ubuntu:24.04

ARG SEQKIT_VERSION=2.13.0
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates wget \
    && wget -qO- "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz" \
    | tar xz -C /usr/local/bin \
    && chmod +x /usr/local/bin/seqkit \
    && apt-get purge -y wget \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser
USER appuser

ENTRYPOINT ["/usr/bin/env"]
CMD ["bash"]
