# BaseJumper BASEJ-SOMATIC

The BioSkryb BASEJ-SOMATIC pipeline performs single-cell somatic SNP/indel
detection and heuristic-QC filtering. Starting from per-sample germline VCFs, it
merges calls across a group of cells, runs a pileup-based heuristic filter and
binomial/beta-binomial statistical testing (Sequoia), removes germline/bulk
variants, optionally annotates with VEP, and produces per-group variant tables,
provenance reports, VAF/digital heatmaps and a MultiQC summary.

## Starting from a VCF (open-source run — no license required)

This README documents the **open-source** run, which **starts from a pre-computed
per-sample VCF** supplied via the `vcf` column of the input CSV.

The only proprietary step in the pipeline is `SENTIEON_DNASCOPE` (Sentieon
DNAscope germline calling, BAM/CRAM → VCF), which **requires a Sentieon license**.
When you provide the `vcf` column, DNAscope is **skipped entirely** and no license
is needed — the pipeline consumes your VCFs directly. (DNAscope only runs for
samples that have no `vcf` value; that path is proprietary and out of scope for
the public release.)

The alignment file (BAM/CRAM) is still required alongside the VCF because the
downstream pileup/heuristic filter reads it.

All custom container images for the open-source run are **built locally** from the
Dockerfiles in [`container/`](container/README.md). Public biocontainers
(bcftools, samtools, ensembl-vep, multiqc) are pulled directly from `quay.io`.

# Pipeline Overview

Main stages of the open-source (VCF-start) run:

- Pre-process and normalise each per-sample VCF (**bcftools** / variantannotation)
- Merge per-sample VCFs into a group VCF (**bcftools**)
- Per-chromosome pileup of the chosen positions from the BAM/CRAM (**samtools**)
- Build NV/NR (variant/reference read) matrices and run the Sequoia
  binomial/beta-binomial statistical filter (custom R)
- Heuristic artifact QC (alignment score, clipped-read fraction, base-position,
  strand support) and second-pass Sequoia filtering (custom R)
- Remove germline variants (statistical + optional bulk VCF + optional VEP AF filter)
- ADO (allele drop-out) analysis (custom R)
- (Optional) VEP annotation when `--vep_cache_dir` is set (**ensembl-vep**)
- Per-group variant tables, provenance, VAF-split diagnostics and a **MultiQC** report

# Running Locally

Instructions for running BASEJ-SOMATIC on a local Ubuntu server.

## Install Java

```
sudo apt-get install default-jdk
java -version
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

## Install Docker

See the Docker installation docs for your distribution
(https://docs.docker.com/engine/install/ubuntu/).

## Build the container images

The custom images are built locally (x86/amd64 only) from [`container/`](container/README.md):

```
bash container/build_all_x86.sh
```

> **No Sentieon license required for the VCF-start run.** The proprietary
> `SENTIEON_DNASCOPE` image is not built and is never invoked when the `vcf`
> column is supplied.

# Reference Data

The pipeline needs a reference-genome bundle (FASTA, index). These are read from
the location given by `--genomes_base`, which defaults to the BioSkryb shared S3
path (`s3://bioskryb-shared-data`). For a local run, download the bundle once and
point `--genomes_base` at your local copy:

```
--genomes_base /path/to/local/genomes
```

VEP annotation (optional) additionally needs a VEP cache, passed via
`--vep_cache_dir`. Leave it unset to skip VEP.

# Test Pipeline Execution

**Input CSV (Illumina, VCF-start)** — columns `biosampleName`, `bam`, `group`, `vcf`:

```
biosampleName,bam,group,vcf
cell1,s3://.../cell1.bam,groupA,s3://.../cell1.vcf.gz
cell2,s3://.../cell2.bam,groupA,s3://.../cell2.vcf.gz
```

**Input CSV (Ultima, VCF-start)** — columns `biosampleName`, `cram`, `group`, `vcf`
(optional `crai`):

```
biosampleName,cram,group,vcf
cell1,s3://.../cell1.cram,groupA,s3://.../cell1.vcf.gz
```

**Command (Illumina, open-source VCF-start run)**

```
nextflow run main.nf \
  --platform Illumina \
  --input_csv $PWD/tests/data/inputs/input_somatic.csv \
  --genomes_base /path/to/local/genomes \
  --genome GRCh38 \
  --outputDir test \
  --max_cpus 8 --max_memory 30.GB
```

> The public release is **x86 (amd64) only**.

## Command Options

```
    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

        [required]
        --input_csv     FILE    Path to input CSV file
        --platform      STR     Input platform: Illumina (BAM) or Ultima (CRAM)

        [starting point]
        --input_csv 'vcf' column    Provide a pre-computed per-sample VCF to skip
                                    the proprietary Sentieon DNAscope step (open-source run).

        [optional]
        --genomes_base  STR     Base path to the reference-genome bundle.
                                DEFAULT: s3://bioskryb-shared-data
        --genome        STR     Reference genome. Options: GRCh38, GRCm39
                                DEFAULT: GRCh38
        --gender        STR     male | female (chrX/Y handling). DEFAULT: male
        --bulk_vcf      FILE    Optional bulk VCF for early germline removal. DEFAULT: none
        --vep_cache_dir STR     VEP cache dir; unset skips VEP annotation.
        --outputDir     DIR     Path to run output directory. DEFAULT: results
        --help          BOOL    Display help message
```

## Tool versions (open-source run)

- `bcftools: 1.21` and `1.14`
- `samtools: 1.21`
- `ensembl-vep: 115.2` (annotation) / `111` (variantannotation image)
- `MultiQC: 1.33`
- Custom R filtering stack: `custom_snp_somatic_filter_sequoia_feb2026`

## Outputs

Outputs are written to `--outputDir`, including per-group merged/annotated VCFs,
NV/NR matrices, variant filter-provenance tables, ADO summaries, VAF/digital
heatmaps and a MultiQC report.

## Lineage handoff manifest

For the downstream `basej-lineage` pipeline, the run emits (via
`EMIT_LINEAGE_INPUTS_MANIFEST`):

- `index/lineage_inputs.csv` — `group,param,path` rows (multi-group safe) with the
  `nr_matrix`, `nv_matrix`, `binary_matrix`, `mandatory_variants_qc_status` and
  per-group `input_csv` paths.
- `index/lineage_vcfs_<group>.csv` — per-group `biosampleName,vcf,vcf_index` list.

The `nr_matrix`/`nv_matrix`/`binary_matrix` paths point at the
`downstream_inputs/downstream_inputs_<group>/` bundle's **`PlusMandatoryNonEmpty`**
scheme (`<group>_NR_...PlusMandatoryNonEmpty.tsv`, etc.). This is deliberate: the
`ForPhylogeny` scheme is placed in full by the lineage phylogeny step, which leaves
variant placement empty and skips the VAF/digital heatmaps. The fuller
`PlusMandatoryNonEmpty` matrix keeps variants available for placement so the heatmaps
render.

# Testing

## Test Data Access

Test data is stored on Wasabi-backed S3 at `s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/`.

To access the test data:

**Step 1 — Get your access keys**

Retrieve your AWS credentials from BioSkryb support (contact basejumper-support for the access link).

**Step 2 — Set environment variables**

```bash
export AWS_ACCESS_KEY_ID=<provided_access_key>
export AWS_SECRET_ACCESS_KEY=<provided_secret_key>
export AWS_DEFAULT_REGION=us-east-1
```

## Running a Test

Run the pipeline with the provided test input CSV (5 samples, chr22-only, from VCF, ~30 min):

```bash
nextflow run main.nf \
  --input_csv tests/data/inputs/nftest_input.csv \
  --max_cpus 8 --max_memory 24.GB --architecture x86 \
  --genome GRCh38 --chrs "chr22" \
  --outputDir results_test
```

The test input uses pre-computed DNAscope VCFs, so DNAscope calling is skipped. BAMs are still required for the pileup step.

## nf-test (Automated Testing)

Install nf-test (requires Java 11+):

```bash
curl -fsSL https://code.askimed.com/install/nf-test | bash
mv nf-test /usr/local/bin/
```

Run the automated tests:

```bash
# Run all tests
nf-test test

# Run only the GRCh38 test
nf-test test tests/main.nf.test --tag GRCh38
```


# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the methodologies described.
