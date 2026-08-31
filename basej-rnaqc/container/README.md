# BASEJ-RNAQC Containers (Public Release)

This directory contains the Docker build definitions for every custom image used
by the `basej-rnaqc` pipeline. These images are **not hosted in any registry** —
you build them **locally** from the Dockerfiles here, then reference the resulting
local image names in `nextflow.config`.

Each Dockerfile builds to a local tag of the form `basejumper_<tag>` (e.g.
`basejumper_custom_r_nf_rnaseq_0.10`).

## Images by pipeline stage

| Container folder | Image tag | Process(es) | Notes |
|---|---|---|---|
| `ubuntu/` | `basejumper_ubuntu_24.04_stable` | MERGE_MULTILANE_FASTQ | Minimal ubuntu for cat-based merge |
| `seqkit/` | `basejumper_seqkit-2.13.0` | SEQKIT_SAMPLE | Read subsampling |
| `custom_r_nf_rnaseq/` | `basejumper_custom_r_nf_rnaseq_0.10` | CREATE_QC_REPORT, CREATE_HTSEQ_SUMMARY, MERGE_HTSEQ_SUMMARY, PLOTTER_PCAHEATMAP | R + Bioconductor for gene annotation/parsing and PCA/heatmap |
| `custom_r_qcplots/` | `basejumper_custom_r_qcplots_0.2.3` | RNA_QC_PLOTS, CREATE_HTSEQ_MATRIX | R + Python for QC scoring, Parquet, composition plots |
| `gene_body_coverage/` | `basejumper_gene_body_coverage_0.2.1` | GENE_BODY_COVERAGE_RNA, GENE_BODY_COVERAGE_RNA_PLOT | Python + RSeQC for gene body coverage analysis |

## Public biocontainers (no Dockerfile needed)

These stages use unmodified, publicly available biocontainers pulled directly
from `quay.io`:

| Process | Image |
|---|---|
| SAMTOOLS_SUBSAMPLE_CRAM | `quay.io/biocontainers/samtools:1.21--h50ea8bc_0` |
| FASTP_TRIM | `quay.io/biocontainers/fastp:0.20.1--h8b12597_0` |
| STAR_ALIGN | `quay.io/biocontainers/star:2.7.6a--0` |
| SAMTOOLS_INDEX_FILTER | `quay.io/biocontainers/samtools:1.21--h50ea8bc_0` |
| QUALIMAP_BAMRNA | `quay.io/biocontainers/qualimap:2.2.2d--1` |
| HTSEQ_COUNTS | `quay.io/biocontainers/htseq:0.13.5--py38h803c66d_1` |
| MULTIQC | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |

## Architecture notes

- The public release is **x86 (amd64) only**. All images build for x86; there are
  no arm variants.

## Building an image

Build each image locally from its folder. Each Dockerfile header lists its exact
build command. General form (local tags, no registry):

```bash
docker build -t basejumper_<tag> container/<folder>/
```

For example:

```bash
docker build -t basejumper_custom_r_nf_rnaseq_0.10 container/custom_r_nf_rnaseq/
docker build -t basejumper_custom_r_qcplots_0.2.3 container/custom_r_qcplots/
docker build -t basejumper_gene_body_coverage_0.2.1 container/gene_body_coverage/
```

To build all images at once:

```bash
bash container/build_all_x86.sh
```

After building, point the `container = "..."` entries in `nextflow.config` at the
local image names you built (see `PUBLIC_RELEASE_INSTRUCTIONS.md`).
