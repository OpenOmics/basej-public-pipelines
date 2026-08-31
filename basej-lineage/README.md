# BaseJumper BASEJ-LINEAGE

The BioSkryb BASEJ-LINEAGE pipeline reconstructs single-cell somatic **lineage
(phylogeny)** and **mutational signatures** from pre-built variant matrices. It is
a matrix-consumer pipeline: given paired NR/NV (variant/reference read) matrices,
a binarized variant×sample matrix, and per-sample annotated VCFs, it builds
phylogenetic trees (Sequoia/treemut), places variants on the tree, draws VAF /
digital heatmaps, and derives COSMIC mutational-signature activities
(SigProfiler, plus optional MuSiCaL and SigDyn).

This pipeline uses **only open-source tools** and public/locally-built container
images — **no license is required**. It is typically run downstream of
`basej-somatic`, which emits the NR/NV/binary matrices and annotated VCFs it
consumes.

All custom container images are **built locally** from the Dockerfiles in
[`container/`](container/README.md). Public biocontainers (bcftools, multiqc) are
pulled directly from `quay.io`.

# Pipeline Overview

Two branches run from the input matrices:

**Lineage / phylogeny**
- Build per-sample genotype tables from annotated VCFs (**bcftools**)
- Build phylogenetic trees for SNV / indel / both from NR/NV matrices (Sequoia + treemut, custom R)
- Place variants on the tree and draw VAF + digital heatmaps (custom R + poppler)
- Gather lineage figures into a report bundle

**Mutational signatures**
- Export a somatic variant table from the binary matrix
- Generate SBS/DBS/ID mutation matrices (**SigProfilerMatrixGenerator**)
- Assign COSMIC signature activities (**SigProfilerAssignment**), merge + plot
- Optional: MuSiCaL (Park lab) and SigDyn (Goncalves lab) refitting
- Aggregate into a **MultiQC** report

# Running Locally

## Install Java, AWS CLI, Nextflow, Docker

```
sudo apt-get install default-jdk
```
```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip && sudo ./aws/install
```
```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```
Docker: see https://docs.docker.com/engine/install/ubuntu/.

## Build the container images

The custom images are built locally (x86/amd64 only) from [`container/`](container/README.md):

```
bash container/build_all_x86.sh
```

> **No license required.** This pipeline uses no proprietary tools.

# Inputs

## Using the lineage manifest from basej-somatic

When `basej-somatic` completes, it produces a manifest CSV at
`index/lineage_inputs.csv` containing all the paths needed by this pipeline. You
can pass that file directly via `--lineage_csv` and the pipeline will resolve
`nr_matrix`, `nv_matrix`, `binary_matrix`, `mandatory_variants_qc_status`, and
`input_csv` from it automatically:

```
nextflow run main.nf --lineage_csv /path/to/lineage_inputs.csv
```

The manifest format is `group,param,path` (one row-set per group). The pipeline
uses the first group and warns if additional groups are present.

**Precedence:** Explicit CLI params always override manifest values. For example,
passing `--lineage_csv ... --nr_matrix /my/override.tsv` will use the override for
`nr_matrix` while still resolving the other params from the manifest.

## Direct params (standalone invocation)

Alternatively, provide the pre-built matrices and per-sample VCFs directly (these are
produced by the `basej-somatic` pipeline):

- `--nr_matrix` / `--nv_matrix` — paired NR/NV matrices (rows = `CHROM_POS_REF_ALT`, cols = samples). **Required.** Use the `basej-somatic` **`PlusMandatoryNonEmpty`** scheme from the `downstream_inputs/downstream_inputs_<group>/` bundle (`<group>_NR_...PlusMandatoryNonEmpty.tsv`), **not** the `ForPhylogeny` scheme: the pre-selected `ForPhylogeny` matrix is placed in full by the phylogeny step, which leaves the variant-placement step empty and skips the VAF/digital heatmaps.
- `--binary_matrix` — binarized variant×sample matrix (0/1), same `PlusMandatoryNonEmpty` scheme. Required for the mutational-signature branch (skipped when unset).
- Per-sample annotated VCFs, via **either**:
  - `--input_csv` — CSV with header `biosampleName,vcf,vcf_index`, or
  - `--vcf_dir` — a directory of `*_somatic_annotated.vcf.gz` (+ `.tbi`) pairs.
- `--genotype_bin` (optional) — pre-binarized genotype matrix (0/0.5/1); when set, VAF discretization is skipped.

# Reference Data

Reference genome bundle is read from `--genomes_base` (default
`s3://bioskryb-shared-data`). For a local run, point it at your local copy:

```
--genomes_base /path/to/local/genomes
```

# Test Pipeline Execution

## Via lineage_csv (using basej-somatic output)

```
nextflow run main.nf \
  --lineage_csv /path/to/somatic_output/index/lineage_inputs.csv \
  --cohort_id   cohort \
  --gender      female \
  --genome      GRCh38 \
  --outputDir   test
```

## Via direct params

```
nextflow run main.nf \
  --nr_matrix   $PWD/tests/data/inputs/NR_matrix.tsv \
  --nv_matrix   $PWD/tests/data/inputs/NV_matrix.tsv \
  --binary_matrix $PWD/tests/data/inputs/binary_matrix.tsv \
  --input_csv   $PWD/tests/data/inputs/input_lineage.csv \
  --cohort_id   cohort \
  --gender      female \
  --genome      GRCh38 \
  --outputDir   test \
  --max_cpus 8 --max_memory 30.GB
```

> The public release is **x86 (amd64) only**.

## Command Options

```
    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

        [auto-resolve]
        --lineage_csv   FILE    Manifest CSV from basej-somatic (group,param,path).
                                When provided, resolves nr_matrix/nv_matrix/binary_matrix/
                                mandatory_variants_qc_status/input_csv from the manifest.

        [required — unless resolved via --lineage_csv]
        --nr_matrix     FILE    Pre-built NR (reference read) matrix TSV
        --nv_matrix     FILE    Paired NV (variant read) matrix TSV
        --input_csv     FILE    CSV (biosampleName,vcf,vcf_index)   [or use --vcf_dir]

        [optional]
        --binary_matrix FILE    Binarized variant x sample matrix (enables mutsig branch)
        --vcf_dir       DIR     Directory of *_somatic_annotated.vcf.gz (+ .tbi)
        --genotype_bin  FILE    Pre-binarized genotype matrix (0/0.5/1)
        --mandatory_variants_qc_status FILE  QC status TSV (failing variants removed from heatmaps)
        --cohort_id     STR     Cohort label. DEFAULT: cohort
        --gender        STR     male | female. DEFAULT: female
        --sequoia_phylogeny_mode  STR   snv,indel,both (and/or all). DEFAULT: snv,indel,both
        --musical_enabled   BOOL    Run MuSiCaL refitting. DEFAULT: true
        --sigdyn_enabled    BOOL    Run SigDyn refitting. DEFAULT: true
        --musicatk_enabled  BOOL    Run musicatk (needs R>=4.4; not in default image). DEFAULT: false
        --genome        STR     Reference genome. DEFAULT: GRCh38
        --outputDir     DIR     Path to run output directory. DEFAULT: results
        --help          BOOL    Display help message
```

## Tool versions (open-source run)

- `bcftools: 1.14`
- `MultiQC: 1.33`
- Phylogeny / signatures R+py stack: `custom_snp_somatic_filter_sequoia_feb2026`
- SigProfiler: `sigprofiler-0.1.0`

## Outputs

Outputs are written to `--outputDir`, including phylogenetic trees, variant
placements, VAF/digital heatmaps, COSMIC signature-activity tables and plots, and
a MultiQC report.

# Testing

## Test Data Access

Test data for this pipeline (NR/NV matrices, binary matrix, VCFs) is included
directly in `tests/data/inputs/` as these files are small. The VCFs and matrices
are derived from a basej-somatic run on chr22-only data.

For pipelines that reference S3-hosted test data:

**Step 1 — Get your access keys**

Retrieve your AWS credentials from BioSkryb support (contact basejumper-support for the access link).

**Step 2 — Set environment variables**

```bash
export AWS_ACCESS_KEY_ID=<provided_access_key>
export AWS_SECRET_ACCESS_KEY=<provided_secret_key>
export AWS_DEFAULT_REGION=us-east-1
```

## Running a Test

Run the pipeline with the provided test inputs:

```bash
nextflow run main.nf \
  --vcf_dir tests/data/inputs/vcfs \
  --nr_matrix tests/data/inputs/nftest_NR_matrix.tsv \
  --nv_matrix tests/data/inputs/nftest_NV_matrix.tsv \
  --binary_matrix tests/data/inputs/nftest_binary_matrix.tsv \
  --max_cpus 8 --max_memory 24.GB \
  --genome GRCh38 --gender male --cohort_id test_cohort \
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

# Run only the GRCh38 test
nf-test test tests/main.nf.test --tag GRCh38
```


# Need Help?

If you need any help, please [submit a helpdesk ticket](https://bioskryb.atlassian.net/servicedesk/customer/portal/3/group/14/create/156).

# References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the methodologies described.
