# BaseJumper BASEJ-RNA-QC

The BioSkryb BASEJ-RNA-QC pipeline evaluates the quality of single-cell RNA
libraries and provides QC metrics to assess sequencing quality, gene expression,
and genomic composition.

Users can sequence single-cell RNA libraries and use this pipeline to assess
library quality via STAR alignment, HTSeq gene quantification, Qualimap genomic
proportion analysis, and composite QC scoring.

# Pipeline Overview

Steps and tools used:

- Subsample reads to 2 million using **SeqKit** (FASTQ) to compare metrics across samples
- Evaluate sequencing quality and trim/clip reads using **fastp** (Illumina FASTQ)
- Map reads to the reference genome using **STAR** (2-pass mode)
- Filter to primary alignments and index using **samtools**
- Quantify gene expression using **HTSeq** (gene-level counting)
- Analyze genomic origin of reads (exonic/intronic/intergenic) using **Qualimap**
- Generate PCA and heatmap visualizations of gene expression
- Compute gene body coverage profiles using **RSeQC**
- Generate per-biosample QC plots, consensus scores, and Parquet outputs (custom R/Python)
- Aggregate metrics across biosamples into an overall report using **MultiQC**

The custom container images are **built locally** from the Dockerfiles in
[`container/`](container/README.md) (see that README for the full image list and
build instructions). Public biocontainers (fastp, STAR, samtools, qualimap,
htseq, multiqc) are pulled directly from `quay.io`.

# Running Locally

Instructions for running BASEJ-RNA-QC on a local Ubuntu server.

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

## Resources Required

For a typical RNA-QC dataset, 16 CPU cores and 60 GB of memory are sufficient for
the STAR alignment stage (loads the genome index into memory). You can constrain
resources with:

```
--max_cpus 16 --max_memory 60.GB
```

# Reference Data

The pipeline needs a reference-genome bundle (STAR index, GTF annotation,
transcript-to-gene mapping, and housekeeping gene BED). These are read from the
location given by `--genomes_base`, which defaults to the BioSkryb shared S3 path
(`s3://bioskryb-shared-data`).

For a public/local run, download the reference bundle for your genome once and
point `--genomes_base` at your local copy. The directory you pass must contain the
expected `genomes/...` layout. Then pass it on every run:

```
--genomes_base /path/to/local/genomes
```

# Test Pipeline Execution

**Command (Illumina FASTQ)**

```
nextflow run main.nf \
  --input_csv $PWD/tests/data/inputs/nftest_input.csv \
  --genomes_base /path/to/local/genomes \
  --genome GRCh38 \
  --outputDir test \
  --max_cpus 16 --max_memory 60.GB
```

## Input Options

The input is passed via an `input.csv` file. The platform is auto-detected from
the columns present.

- **Illumina / Element (FASTQ)** — columns `biosampleName`, `read1`, `read2`:

```
biosampleName,read1,read2,groups
sample1,s3://.../sample1_R1.fastq.gz,s3://.../sample1_R2.fastq.gz,Group1
sample2,s3://.../sample2_R1.fastq.gz,s3://.../sample2_R2.fastq.gz,Group1
```

  Multi-lane inputs are supported by pipe-delimiting (`|`) multiple paths within
  `read1`/`read2`.

**Optional `groups` column**: carries sample group information used by the QC
plotting / consensus-score grouping.

## Command Options

```
    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

        [required]
        --input_csv         FILE    Path to input csv file

        [optional]
        --genomes_base      STR     Base path to the reference-genome bundle.
                                    DEFAULT: s3://bioskryb-shared-data

        --genome            STR     Reference genome. Options: GRCh38, GRCm39
                                    DEFAULT: GRCh38

        --outputDir         DIR     Path to run output directory
                                    DEFAULT: results

        --n_reads           VAL     Number of reads to subsample for QC
                                    DEFAULT: 2000000

        --help              BOOL    Display help message
```

## Tool versions

- `SeqKit: 2.13.0`
- `fastp: 0.20.1`
- `STAR: 2.7.6a`
- `samtools: 1.21`
- `Qualimap: 2.2.2d`
- `HTSeq: 0.13.5`
- `RSeQC: 4.0.0`
- `MultiQC: 1.33`

## Outputs

Outputs are written to the directory specified by `--outputDir`, including
per-biosample QC Parquet files, merged metric tables, gene count matrices,
QC plots, consensus scores, and a MultiQC report.


# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the
methodologies described.
