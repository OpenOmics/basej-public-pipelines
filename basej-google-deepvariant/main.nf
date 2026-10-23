nextflow.enable.dsl=2

// ============================================================================
// BASEJ-GOOGLE-DEEPVARIANT: Google DeepVariant Variant Calling Pipeline
// ============================================================================
// Description: Runs after basej-wgs. Takes per-biosample BAM or CRAM input
//              and produces small-variant VCFs using the three-stage Google
//              DeepVariant flow:
//                make_examples (CPU, sharded)
//                call_variants (GPU)
//                postprocess_variants (CPU)
//              DeepVariant natively supports CRAM via htslib (--ref decodes
//              CRAM automatically since DV v0.9.0), so no CRAM→BAM conversion
//              is needed.
//              Uses BioSkryb's custom DV model + population VCFs from
//              genomes.config.
// Outputs: Per-biosample <sample>_deepvariant.vcf.gz{,.tbi}
//          + bcftools stats + MultiQC summary report
// ============================================================================
//
// Usage:
//   nextflow run main.nf \
//     --input_csv samples.csv \
//     --genome GRCh38 \
//     --architecture arm \
//     -profile batch_dev
//
// Input CSV columns (one of bam | cram per row):
//   biosampleName, bam[, bai]
//   biosampleName, cram[, crai]
// ============================================================================


// ============================================================================
// PROCESS: DEEPVARIANT_MAKE_EXAMPLES_ONLY  (inline)
// Description: Stage 1 of DeepVariant. CPU-bound; shards example generation
//              across task.cpus using GNU parallel. Accepts BAM or CRAM input
//              natively (DeepVariant uses htslib and decodes CRAM via --ref).
// ============================================================================
process DEEPVARIANT_MAKE_EXAMPLES_ONLY {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(alignment), path(alignment_index)
    path(deepvariant_model)
    path(population_vcfs)
    path(reference)
    path(regions)

    output:
    tuple val(sample_name), path("${sample_name}.tfrecord-*-of-*.gz"), val(task.cpus), emit: example_tfrecords
    // gVCF (non-variant site) tfrecords — only emitted when --make_gvcf is set.
    // Sharded across the same task.cpus as the example tfrecords; the shard
    // count (task.cpus) is carried through so postprocess can reconstruct the
    // --nonvariant_site_tfrecord_path @<shards> pattern.
    tuple val(sample_name), path("${sample_name}.gvcf.tfrecord-*-of-*.gz"), val(task.cpus), emit: gvcf_tfrecords, optional: true

    script:
    def regions_arg = params.mode == "exome" ? "--regions ${regions}" : ""
    def gvcf_arg    = params.make_gvcf ? "--gvcf \"${sample_name}.gvcf.tfrecord@${task.cpus}.gz\"" : ""
    """
    set -e
    export DV_BIN_PATH=/opt/deepvariant/bin
    mkdir -p scratch_${sample_name}
    export TMPDIR=\$PWD/scratch_${sample_name}

    population_vcfs_list=\$(ls ${population_vcfs}/cohort-chr*.release_missing2ref.no_calls.vcf.gz 2>/dev/null | tr '\\n' ',' | sed 's/,\$//' || echo "")

    seq 0 ${task.cpus - 1} | \\
    parallel -q --halt 2 --line-buffer \\
    make_examples \\
      --mode calling \\
      --ref ${reference}/genome.fa \\
      --reads ${alignment} \\
      --examples "${sample_name}.tfrecord@${task.cpus}.gz" \\
      --checkpoint "${deepvariant_model}/${params.checkpoint_filename}" \\
      --population_vcfs="\${population_vcfs_list}" \\
      ${gvcf_arg} \\
      ${regions_arg} \\
      --task {}
    """
}


// ============================================================================
// PROCESS: DEEPVARIANT_CALL_VARIANTS  (inline — GPU)
// Description: Stage 2 of DeepVariant. Runs the trained model on the example
//              tfrecords. GPU-accelerated container.
// ============================================================================
process DEEPVARIANT_CALL_VARIANTS {
    label "gpu"
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(example_tfrecords), val(shards)
    path(deepvariant_model)

    output:
    tuple val(sample_name), path("${sample_name}_variants_output*tfrecord.gz"), path("variants_shards.txt"), emit: variants_output

    script:
    """
    set -e

    call_variants \\
      --examples=${sample_name}.tfrecord@${shards}.gz \\
      --outfile=${sample_name}_variants_output.tfrecord.gz \\
      --checkpoint="${deepvariant_model}/${params.checkpoint_filename}" \\
      --batch_size=${params.call_variants_batch_size}

    # Count files produced (used by postprocess_variants' shard pattern)
    ls -1 ${sample_name}_variants_output*.tfrecord.gz | wc -l > variants_shards.txt
    """
}


// ============================================================================
// PROCESS: DEEPVARIANT_POSTPROCESS  (inline)
// Description: Stage 3 of DeepVariant. Assembles the final VCF + index from
//              the per-shard call_variants output.
// ============================================================================
process DEEPVARIANT_POSTPROCESS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(variants_output), path(variants_shards_file), path(gvcf_tfrecords), val(gvcf_shards)
    path(reference)

    output:
    tuple val(sample_name), path("${sample_name}_deepvariant.vcf.gz"), path("${sample_name}_deepvariant.vcf.gz.tbi"), emit: vcf
    // gVCF (+ index) — only emitted when --make_gvcf is set.
    tuple val(sample_name), path("${sample_name}_deepvariant.g.vcf.gz"), path("${sample_name}_deepvariant.g.vcf.gz.tbi"), emit: gvcf, optional: true
    path("deepvariant_version.yml"),                                                                                  emit: version

    script:
    def gvcf_arg = params.make_gvcf ? "--gvcf_outfile=${sample_name}_deepvariant.g.vcf.gz --nonvariant_site_tfrecord_path=${sample_name}.gvcf.tfrecord@${gvcf_shards}.gz" : ""
    """
    set -e

    variants_output_shards=\$(cat ${variants_shards_file} | tr -d '[:space:]')
    if [ -z "\$variants_output_shards" ] || [ "\$variants_output_shards" -lt 1 ]; then
        echo "Error: invalid shard count: \$variants_output_shards" >&2
        exit 1
    fi

    postprocess_variants \\
      --ref=${reference}/genome.fa \\
      --infile=${sample_name}_variants_output@\${variants_output_shards}.tfrecord.gz \\
      --outfile=${sample_name}_deepvariant.vcf.gz \\
      ${gvcf_arg} \\
      --cpus=${task.cpus}

    # Make sure .tbi exists; postprocess writes it but we tabix as a fallback.
    if [ ! -s "${sample_name}_deepvariant.vcf.gz.tbi" ]; then
        echo "tbi missing, regenerating with tabix"
        tabix -p vcf ${sample_name}_deepvariant.vcf.gz
    fi

    # Same fallback for the gVCF index when gVCF output is enabled.
    if [ "${params.make_gvcf}" = "true" ] && [ -s "${sample_name}_deepvariant.g.vcf.gz" ] && [ ! -s "${sample_name}_deepvariant.g.vcf.gz.tbi" ]; then
        echo "gvcf tbi missing, regenerating with tabix"
        tabix -p vcf ${sample_name}_deepvariant.g.vcf.gz
    fi

    echo "DeepVariant: ${params.deepvariant_version}" > deepvariant_version.yml
    """
}


// ============================================================================
// PROCESS: BCFTOOLS_STATS  (inline)
// Description: Per-sample VCF stats for MultiQC. Lightweight wrapper around
//              `bcftools stats` so the report can show SNV/indel counts,
//              ts/tv, etc., per biosample.
// ============================================================================
process BCFTOOLS_STATS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(vcf), path(tbi)

    output:
    path("${sample_name}.bcftools.stats.txt"), emit: stats

    script:
    """
    set -e
    bcftools stats --threads ${task.cpus} ${vcf} > ${sample_name}.bcftools.stats.txt
    """
}


// ============================================================================
// PROCESS: MULTIQC_DV  (inline)
// Description: Aggregate per-sample bcftools stats + DeepVariant version
//              files into a single MultiQC HTML report.
// ============================================================================
process MULTIQC_DV {
    input:
    path(input_files)
    val(project)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)

    output:
    path("multiqc_report.html"),  emit: report
    path("multiqc_report_data"),  emit: data

    script:
    """
    cat > multiqc_config.yaml << EOF
custom_logo_title: 'BioSkryb Genomics'
title: "basej-google-deepvariant v${pipeline_version}"
report_header_info:
  - Project: "${project}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false

run_modules:
  - bcftools
  - custom_content
EOF

    multiqc . \\
        --config multiqc_config.yaml \\
        --filename multiqc_report \\
        --force
    """
}


// ============================================================================
// MAIN WORKFLOW
// ============================================================================
workflow {
    main:

    if (!params.input_csv) {
        exit 1, "ERROR: --input_csv parameter is required"
    }

    if (params.genome != "GRCh38") {
        log.warn "basej-google-deepvariant has only been validated on GRCh38; --genome=${params.genome} may not have a matching DV model or population VCFs."
    }

    // Exome mode requires a regions BED — the schema enforces this for the
    // launcher (Tower / MCP), but a CLI run with --mode exome and no
    // --regions would otherwise silently call genome-wide. Fail fast.
    if (params.mode == 'exome' && (!params.regions || params.regions.toString().trim() == "")) {
        exit 1, "ERROR: --regions is required when --mode=exome (BED of calling intervals). Set --regions /path/to/exome.bed or switch to --mode=wgs."
    }

    // -------------------------------------------------------------------------
    // Phase 0: Parse input CSV
    //   biosampleName, bam[, bai]   — BAM input rows
    //   biosampleName, cram[, crai] — CRAM input rows (passed directly to DV)
    // -------------------------------------------------------------------------
    ch_input = Channel.fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def name = row.biosampleName ?: row.sample_name ?: row.sample
            if (!name) { exit 1, "ERROR: input CSV row is missing biosampleName: ${row}" }

            if (row.cram) {
                def cram = file(row.cram, checkIfExists: true)
                def crai = row.crai ? file(row.crai, checkIfExists: true) : file("${row.cram}.crai", checkIfExists: true)
                tuple(name, cram, crai)
            } else if (row.bam) {
                def bam = file(row.bam, checkIfExists: true)
                def bai = row.bai ? file(row.bai, checkIfExists: true) : file("${row.bam}.bai", checkIfExists: true)
                tuple(name, bam, bai)
            } else {
                exit 1, "ERROR: input CSV row must contain either 'bam' or 'cram': ${row}"
            }
        }

    ch_input.ifEmpty { exit 1, "ERROR: No samples found in --input_csv" }

    // -------------------------------------------------------------------------
    // Phase 1: Three-stage DeepVariant (BAM and CRAM accepted natively)
    // -------------------------------------------------------------------------
    ch_deepvariant_model = file(params.deepvariant_model, checkIfExists: true)
    ch_population_vcfs   = file(params.population_vcfs,   checkIfExists: true)
    ch_reference         = file(params.reference,         checkIfExists: true)
    ch_regions           = params.regions ? file(params.regions, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")

    DEEPVARIANT_MAKE_EXAMPLES_ONLY(
        ch_input,
        ch_deepvariant_model,
        ch_population_vcfs,
        ch_reference,
        ch_regions
    )

    DEEPVARIANT_CALL_VARIANTS(
        DEEPVARIANT_MAKE_EXAMPLES_ONLY.out.example_tfrecords,
        ch_deepvariant_model
    )

    // Build the postprocess input channel.
    //   call_variants output : [ name, variants_output, variants_shards_file ]
    //   make_examples gvcf   : [ name, gvcf_tfrecords, gvcf_shards ]
    // The gVCF tfrecords are produced by make_examples (NOT call_variants) and
    // must be carried to postprocess via a join on sample_name so each sample's
    // gVCF tfrecords stay matched with its own variant calls. When gVCF is
    // disabled we pad with a placeholder file + zero shard count to keep the
    // process input arity uniform.
    if (params.make_gvcf) {
        ch_postprocess_in = DEEPVARIANT_CALL_VARIANTS.out.variants_output
            .join(DEEPVARIANT_MAKE_EXAMPLES_ONLY.out.gvcf_tfrecords, by: 0)
            // -> [ name, variants_output, variants_shards_file, gvcf_tfrecords, gvcf_shards ]
    } else {
        ch_postprocess_in = DEEPVARIANT_CALL_VARIANTS.out.variants_output
            .map { name, variants_output, variants_shards_file ->
                tuple(name, variants_output, variants_shards_file, file("${projectDir}/assets/NO_FILE"), 0)
            }
    }

    DEEPVARIANT_POSTPROCESS(
        ch_postprocess_in,
        ch_reference
    )

    // -------------------------------------------------------------------------
    // Phase 3: Per-sample bcftools stats + MultiQC summary
    // -------------------------------------------------------------------------
    BCFTOOLS_STATS(DEEPVARIANT_POSTPROCESS.out.vcf)

    ch_multiqc_inputs = BCFTOOLS_STATS.out.stats
        .mix(DEEPVARIANT_POSTPROCESS.out.version.first())
        .collect()

    MULTIQC_DV(
        ch_multiqc_inputs,
        params.project,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version ?: '1.0.0'
    )

    // -------------------------------------------------------------------------
    // Publish targets (see output{} block below). Without this section the
    // pipeline produced results only inside the work dir and nothing landed in
    // workflow_outputs/, so resumed runs appeared to have "no output".
    // -------------------------------------------------------------------------
    publish:
    // Per-biosample DeepVariant VCFs (+ tbi). Indexed so downstream steps can
    // locate each biosample's VCF by name.
    vcf_files = DEEPVARIANT_POSTPROCESS.out.vcf
        .map { sample_name, vcf, tbi -> [biosampleName: sample_name, vcf: vcf, tbi: tbi] }

    // Per-biosample DeepVariant gVCFs (+ tbi). Only populated when --make_gvcf
    // is set; otherwise this channel is empty and nothing is published.
    gvcf_files = DEEPVARIANT_POSTPROCESS.out.gvcf
        .map { sample_name, gvcf, tbi -> [biosampleName: sample_name, gvcf: gvcf, tbi: tbi] }

    // Per-biosample bcftools stats (feed MultiQC; published for provenance).
    bcftools_stats = BCFTOOLS_STATS.out.stats

    // MultiQC summary report + its data directory.
    multiqc_report = MULTIQC_DV.out.report
    multiqc_data   = MULTIQC_DV.out.data
}


// ============================================================================
// OUTPUT CONFIGURATION
// ============================================================================
// Resolved against the top-level `outputDir` config setting (= params.outputDir).
// Paths mirror the basej-wgs convention so resumed runs republish cached
// outputs to workflow_outputs/<workspace>/<workflow_id>/ for the original id.
output {
    vcf_files {
        path "vcf/${params.workspace}/dna/tool=deepvariant-${params.deepvariant_version}"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/vcf.csv"
            header true
        }
        tags workspace:     params.workspace,
             workflow_id:   params.workflow_id,
             project:       params.project,
             pipeline:      "basej-google-deepvariant",
             molecule_type: "dna",
             artifact:      "vcf",
             tool:          "deepvariant",
             reference:     params.genome
    }

    gvcf_files {
        path "vcf/${params.workspace}/dna/tool=deepvariant-${params.deepvariant_version}/gvcf"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/gvcf.csv"
            header true
        }
        tags workspace:     params.workspace,
             workflow_id:   params.workflow_id,
             project:       params.project,
             pipeline:      "basej-google-deepvariant",
             molecule_type: "dna",
             artifact:      "gvcf",
             tool:          "deepvariant",
             reference:     params.genome
    }

    bcftools_stats {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/bcftools_stats"
        tags workspace:     params.workspace,
             workflow_id:   params.workflow_id,
             project:       params.project,
             pipeline:      "basej-google-deepvariant",
             molecule_type: "dna",
             artifact:      "bcftools_stats",
             tool:          "bcftools",
             reference:     params.genome
    }

    multiqc_report {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/reports"
        tags workspace:   params.workspace,
             workflow_id: params.workflow_id,
             project:     params.project,
             pipeline:    "basej-google-deepvariant",
             artifact:    "multiqc_report"
    }

    multiqc_data {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/reports"
        tags workspace:   params.workspace,
             workflow_id: params.workflow_id,
             project:     params.project,
             pipeline:    "basej-google-deepvariant",
             artifact:    "multiqc_data"
    }
}
