# BaseJumper Public Pipelines

Public pipelines from BioSkryb Genomics for single-cell genomic analysis. Each
pipeline lives in its own top-level subdirectory and is self-contained, with its
own run instructions.

## Pipelines

| Pipeline | Directory | Description |
|---|---|---|
| basej-dnaqc | [`basej-dnaqc/`](basej-dnaqc/) | Single-cell DNA low-pass QC — subsample, align, dedup, CNV, per-biosample QC metrics and consensus scoring |
| basej-rnaqc | [`basej-rnaqc/`](basej-rnaqc/) | Single-cell RNA library QC — STAR 2-pass alignment, HTSeq gene quantification, Qualimap genomic composition, and composite QC scoring |
| basej-wgs | [`basej-wgs/`](basej-wgs/) | Single-cell WGS/WES QC — align, dedup, and collect whole-genome or hybrid-selection coverage metrics |
| basej-google-deepvariant | [`basej-google-deepvariant/`](basej-google-deepvariant/) | Germline SNV/indel calling from BAM/CRAM with Google DeepVariant using a BioSkryb custom-trained model that corrects PTA artifacts |
| basej-somatic | [`basej-somatic/`](basej-somatic/) | Single-cell somatic SNP/indel detection and heuristic-QC filtering with per-group variant tables |
| basej-lineage | [`basej-lineage/`](basej-lineage/) | Single-cell lineage/phylogeny reconstruction and COSMIC mutational-signature analysis from variant matrices |

## Repository layout

```
basej-public-pipelines/
├── README.md          # this file
└── <pipeline>/        # one directory per pipeline
    ├── README.md      # pipeline-specific run instructions and options
    ├── main.nf        # pipeline entrypoint
    ├── container/     # Dockerfiles for the custom images used by the pipeline
    ├── conf/          # pipeline configuration
    └── tests/         # example inputs and tests
```

## Getting started

Pick the pipeline you need from the table above and follow the `README.md` inside
its directory for setup, inputs, run commands, options, and outputs.

To obtain the reference-genome bundles and pipeline resources needed to run these
pipelines, email [basejumper@bioskryb.com](mailto:basejumper@bioskryb.com).

## Running on Biowulf

The Basejumper pipeline suite is hosted under OpenOmics on the Biowulf cluster and can be utilized by loading the module as shown below.

There are 6 individual pipelines available:
- dnaqc, rnaqc, wgs, deepvariant, somatic, and lineage

There is 1 pipeline chain available as an automated option for running a high-pass workflow:
- wgs -> deepvariant -> somatic -> lineage

```
# Load the Basejumper module
module load openomics/basejumper

# Basejumper Pipeline Suite

Usage: basej-run [--help] \
          --pipeline {dnaqc,rnaqc,wgs,deepvariant,somatic,lineage,wgs2lineage} \
          --input_csv INPUT_CSV [--lineage_csv LINEAGE_CSV] \
          --outputDir OUTPUTDIR \
          [--workDir WORKDIR] [--dry-run] [--force] \
          [--wgs_extra WGS_EXTRA] [--deepvariant_extra DEEPVARIANT_EXTRA] \
          [--somatic_extra SOMATIC_EXTRA] [--lineage_extra LINEAGE_EXTRA]

Description:
  Run any individual or the wgs2lineage BioSkryb BaseJumper pipeline.

Required arguments:
  --pipeline {dnaqc,rnaqc,wgs,deepvariant,somatic,lineage,wgs2lineage}
                                  Select which pipeline to run. 'wgs2lineage' is a convenience alias 
                                  for the chain: wgs -> deepvariant -> somatic -> lineage.
                                
  --input_csv INPUT_CSV
                                  Input CSV. Required for --pipeline dnaqc / rnaqc / wgs / deepvariant / somatic / wgs2lineage.
                                  For wgs2lineage, this must be the wgs input samplesheet (biosampleName,read1,read2,group);
                                  basej-run auto-builds the samplesheets for deepvariant and somatic using pipeline-provided index files.

  --lineage_csv LINEAGE_CSV
                                  Input CSV for --pipeline lineage specifically. Output from basej-somatic.
  --outputDir OUTPUTDIR
                                  Base output directory (default: ./basej-results).

Optional arguments:
  --workdir WORKDIR               Nextflow work directory (default: ./work).
  --dry-run                       Print the nextflow run commands.
  --force                         Rerun a step even if already marked complete.
  -h, --help                      Show this help message and exit.

Optional parameter overrides for --pipeline wgs2lineage (quoted as one string):
  --wgs_extra WGS_EXTRA                     Additional basej-wgs parameter flags

  --deepvariant_extra DEEPVARIANT_EXTRA     Additional basej-deepvariant parameter flags

  --somatic_extra SOMATIC_EXTRA             Additional basej-somatic parameter flags

  --lineage_extra LINEAGE_EXTRA             Additional basej-lineage parameter flags

```

### Example Run Commands
```
  # Run a standalone pipeline (dnaqc, rnaqc, wgs, deepvariant, somatic)
  basej-run --pipeline dnaqc --input_csv dnaqc_input.csv --outputDir results [options]

  # Run the lineage pipeline using somatic's output samplesheet
  basej-run --pipeline lineage --lineage_csv results/somatic/index/lineage_inputs.csv --outputDir results [options]

  # Run the chain wgs2lineage (wgs -> deepvariant -> somatic -> lineage)
  basej-run --pipeline wgs2lineage --input_csv wgs_input.csv --outputDir results

  # Run the full chain wgs2lineage with a pipeline-specific parameter overridden
  basej-run --pipeline wgs2lineage --input_csv wgs_input.csv --outputDir results \
    --somatic_extra "--gender female"
```

## Need help?

If you need any help, please email
[basejumper@bioskryb.com](mailto:basejumper@bioskryb.com).

## References

NOTE: Several studies have utilized BaseJumper pipelines as part of standard
quality control processes implemented through ResolveServices<sup>SM</sup>. While
these pipelines may not be explicitly cited, they are integral to the
methodologies described.
