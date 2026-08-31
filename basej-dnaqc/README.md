# BaseJumper BASEJ-DNA-QC

The BioSkryb BASEJ-DNA-QC pipeline evaluates the quality of a single-cell DNA
library and provides several QC metrics to assess the quality of the sequencing
reads.

One way that users can ensure a single-cell library is uniformly amplified with
low allelic dropout is to first sequence using "low-pass" (low-throughput)
sequencing of around 2M reads per sample. Data from the low-pass run can be used
to estimate genome coverage if the libraries were to be used for high-depth
sequencing, so users can select only quality libraries for high-depth runs.

This README documents the **fully open-source** run of the pipeline, which is now
the **default** (`pipeline_tool = opensource`). It uses only open-source tools and
public container images — no proprietary Sentieon license is required. (To use the
proprietary Sentieon path instead, pass `--pipeline_tool sentieon`.)

# Pipeline Overview

Steps and tools used by the open-source run:

- Subsample reads to 2 million using **SeqKit** (FASTQ) or **samtools** (Ultima CRAM) to compare metrics across samples
- Evaluate sequencing quality and trim/clip reads using **fastp** (Illumina FASTQ)
- Map reads to the reference genome using **BWA-MEM2**
- Sort and index alignments using **samtools**
- Remove duplicate reads using **samtools markdup** (`-r`, matching Sentieon's `Dedup --rmdup`, so duplicates are dropped from the BAM rather than only flagged)
- Collect alignment, GC bias, insert size, mean-quality, quality-yield, and coverage metrics using **Picard** (`CollectMultipleMetrics`) plus `samtools depth`
- Estimate library complexity using **preseq** (`bam2mr` + `gc_extrap`)
- Evaluate copy-number variation (CNV) using a custom **Ginkgo** implementation (**bedtools** + Ginkgo R)
- Generate per-biosample QC plots, consensus scores, and Parquet outputs (custom R/Python)
- Aggregate metrics across biosamples and tools into an overall report using **MultiQC**

The custom container images for the open-source run are **built locally** from the
Dockerfiles in [`container/`](container/README.md) (see that README for the full
image list and build instructions). Public biocontainers (fastp, bedtools,
multiqc) are pulled directly from `quay.io`.

# Running Locally

Instructions for running BASEJ-DNA-QC on a local Ubuntu server.

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

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

> **No license required.** The open-source run does not use Sentieon, so there is
> no license setup step. Build the custom images locally from [`container/`](container/README.md);
> public biocontainers (fastp, bedtools, multiqc) are pulled from `quay.io`.

## Resources Required

For a typical low-pass dataset (less than 8 million reads), 8 CPU cores and 30 GB
of memory are sufficient for the alignment stage. You can constrain resources with:

```
--max_cpus 8 --max_memory 30.GB
```

# Reference Data

The pipeline needs a reference-genome bundle (FASTA, BWA-MEM2 index, intervals,
and the Ginkgo CNV references). These are read from the location given by
`--genomes_base`, which defaults to the BioSkryb shared S3 path
(`s3://bioskryb-shared-data`).

For a public/local run, download the reference bundle for your genome once and
point `--genomes_base` at your local copy. The directory you pass must contain the
expected `genomes/...` layout (e.g.
`<genomes_base>/genomes/Homo_sapiens/NCBI/GRCh38/...`). Then pass it on every run:

```
--genomes_base /path/to/local/genomes
```

This avoids per-run S3 access and lets the pipeline run fully offline once the
data and container images are in place.

# Test Pipeline Execution

The example **input** datasets referenced below are publicly available under
`s3://bioskryb-public-data/pipeline_resources` and are fetched automatically
during the run. The **reference genome** is resolved from `--genomes_base` (see
[Reference Data](#reference-data) above).

**Command (Illumina FASTQ, open-source run)**

```
nextflow run main.nf \
  --input_csv $PWD/tests/data/inputs/input_qcTest3.csv \
  --genomes_base /path/to/local/genomes \
  --architecture x86 \
  --genome GRCh38 \
  --outputDir test \
  --max_cpus 8 --max_memory 30.GB
```

**Command (Ultima CRAM, open-source run)**

```
nextflow run main.nf \
  --input_csv $PWD/tests/data/inputs/input_ultima_subsampled.csv \
  --genomes_base /path/to/local/genomes \
  --architecture x86 \
  --genome GRCh38 \
  --outputDir test \
  --max_cpus 8 --max_memory 30.GB
```

> The open-source run is the default, so `--pipeline_tool` is not required. The
> open-source alignment stages (BWA-MEM2, Picard) are x86-only, so use
> `--architecture x86`. The subsampling and merge stages additionally provide
> arm-native images if you run those on Graviton.

## Input Options

The input is passed via an `input.csv` file. The platform is auto-detected from
the columns present.

- **Illumina / Element (FASTQ)** — columns `biosampleName`, `read1`, `read2`:

```
biosampleName,read1,read2
chr22_testsample1,s3://.../chr22_testsample1_R1.fastq.gz,s3://.../chr22_testsample1_R2.fastq.gz
chr22_testsample2,s3://.../chr22_testsample2_R1.fastq.gz,s3://.../chr22_testsample2_R2.fastq.gz
```

  Multi-lane inputs are supported by pipe-delimiting (`|`) multiple paths within
  `read1`/`read2`.

- **Ultima (CRAM)** — columns `biosampleName`, `cram` (the `.cram.crai` index is
  auto-discovered, or supply an optional `crai` column):

```
biosampleName,cram
HG001-UltimaWGS_Z0022,s3://.../HG001-UltimaWGS_Z0022.cram
```

**Optional `groups` column**: an optional `groups` column can carry sample group
information used by the QC plotting / consensus-score grouping.

## Command Options

```
    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

        [required]
        --input_csv         FILE    Path to input csv file

        [pipeline tool]
        --pipeline_tool     STR     Tool selection: 'opensource' (BWA-MEM2 + samtools markdup -r + Picard)
                                    or 'sentieon' (proprietary).
                                    DEFAULT: opensource

        --architecture      STR     Compute architecture: 'arm' (Graviton) or 'x86'.
                                    Use 'x86' for the open-source run (alignment is x86-only).
                                    DEFAULT: arm

        [optional]
        --genomes_base      STR     Base path to the reference-genome bundle. Set this to your
                                    local download directory for public/offline runs.
                                    DEFAULT: s3://bioskryb-shared-data

        --genome            STR     Reference genome. Options: GRCh38, GRCm39, ARSUCD2
                                    DEFAULT: GRCh38

        --outputDir         DIR     Path to run output directory
                                    DEFAULT: results

        --n_reads           VAL     Number of reads to subsample for QC
                                    DEFAULT: 2000000

        --read_length       VAL     Read length used for Ginkgo reference selection
                                    DEFAULT: 50

        --bin_size          VAL     Ginkgo CNV bin size (500000, 1000000, 2000000)
                                    DEFAULT: 1000000

        --help              BOOL    Display help message
```

## Tool versions (open-source run)

- `SeqKit: 2.13.0`
- `fastp: 0.20.1`
- `BWA-MEM2: 2.2.1`
- `samtools: 1.23.1`
- `Picard: 3.1.x`
- `preseq: 2.0.3` (includes `bam2mr`)
- `bedtools: 2.28.0`
- `Ginkgo: 0.3.1`
- `MultiQC: 1.33`

## Outputs

Outputs are written to the directory specified by `--outputDir`, including
per-biosample QC Parquet files, merged metric tables, CNV summaries, QC plots,
consensus scores, and a MultiQC report.


# Docker User / File Permissions

The custom images run as a non-root user (`appuser`). On hosts where your user ID
differs from the image's, containerized tasks can fail to write to the Nextflow work
directory with `touch: cannot touch '.command.trace': Permission denied` (typically
on the Ginkgo and preseq steps). If you hit this, tell Nextflow to run the task
containers as your own user:

1. Create a file named `local.config` next to `main.nf` with this content:

   ```groovy
   docker {
       runOptions = '-u $(id -u):$(id -g)'
   }
   ```

2. Add `-c local.config` to your `nextflow run` command, for example:

   ```
   nextflow run main.nf \
     -c local.config \
     --input_csv $PWD/tests/data/inputs/input_qcTest3.csv \
     --genomes_base /path/to/local/genomes \
     --architecture x86 --genome GRCh38 \
     --outputDir test --max_cpus 8 --max_memory 30.GB
   ```

`$(id -u)`/`$(id -g)` are filled in by your shell at runtime, so the containers run
as you and can write to the work directory.


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

Run the pipeline with the provided test input CSV:

```bash
nextflow run main.nf \
  --input_csv tests/data/inputs/nftest_input.csv \
  --max_cpus 8 --max_memory 24.GB --architecture x86 \
  --genome GRCh38 --n_reads 500000 \
  --pipeline_tool opensource \
  --outputDir results_test
```

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

# Run only the GRCh38 opensource test
nf-test test tests/main.nf.test --tag GRCh38
```

> The `GRCh38` opensource test is the one exercised automatically by the release CI.


# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

- Chung, C., Yang, X., Hevner, R. F., et al. (2024). Cell-type-resolved mosaicism
  reveals clonal dynamics of the human forebrain. Nature, 629(8011), 384–392.
  [https://doi.org/10.1038/s41586-024-07292-5](https://doi.org/10.1038/s41586-024-07292-5)

- Zhao, Y., Luquette, L. J., Veit, A. D., et al. (2024). High-resolution detection
  of copy number alterations in single cells with HiScanner. BioRxiv.
  [https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full](https://www.biorxiv.org/content/10.1101/2024.04.26.587806v1.full)

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While these pipelines may not be explicitly cited, they are integral to the
methodologies described.
