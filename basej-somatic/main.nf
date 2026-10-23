nextflow.enable.dsl=2

// ============================================================================
// BASEJ-SOMATIC: Single-Cell Somatic Variant Calling Pipeline
// ============================================================================
// Description: BAM/CRAM input → DNAscope germline VCF → somatic SNP/indel
//              filtering with heuristic QC scoring. Produces variant tables,
//              phylogenetic trees, VAF heatmaps, and a MultiQC summary report.
//              Supports both Illumina (BAM) and Ultima (CRAM) input platforms.
// Outputs: Per-group variant tables + provenance + MultiQC report
// ============================================================================
//
// Usage:
//   nextflow run main.nf \
//     --platform Illumina --input_csv illumina_samples.csv ...
//   
//   nextflow run main.nf \
//     --platform Ultima --input_csv ultima_samples.csv ...
//
// Input CSV (Illumina): biosampleName, bam, group [, vcf]
// Input CSV (Ultima):   biosampleName, cram, group [, vcf, crai]
// ============================================================================

// Import all somatic-filtering processes from local modules (no publishDir)
include { PREPROCESS_VCF }                                                         from './modules.nf'
include { MERGE_PROCESSED_VCF }                                                    from './modules.nf'
include { CUSTOM_BAM_GROUP_PILEUP }                                                from './modules.nf'
include { CREATE_TAB_NVNR }                                                        from './modules.nf'
include { SEQUOIA_BINOM_BETABINOM_TAB_NV_NR }                                      from './modules.nf'
include { CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR }                                from './modules.nf'
include { CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR } from './modules.nf'
include { CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES }                                from './modules.nf'
include { CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS }                                      from './modules.nf'
include { BULK_GET_VARIANTS_TO_FILTER }                                            from './modules.nf'
include { CREATE_EMPTY_BULK_VARIANTS }                                             from './modules.nf'
include { FILTER_CHOSEN_VARIANTS_BY_BULK }                                         from './modules.nf'
include { GET_VARIANTS_FROM_MERGED_VCF }                                           from './modules.nf'
include { GET_LIST_POS_FROM_CHOSEN_VARIANTS }                                      from './modules.nf'
include { FILTER_DF_NV_BY_CHOSEN_VARIANTS }                                        from './modules.nf'
include { SEQUOIA_SECOND_FILTER }                                                  from './modules.nf'
include { SUBSET_MERGED_VCF_CHOSEN_VARIANTS; SPLIT_SUBSET_VCF_BY_CHR;
          VEP_ANNOTATE; SORT_INDEX_VEP_VCF; MERGE_VEP_VCF_BY_GROUP }              from './modules.nf'
include { FILTER_VEP_GERMLINE }                                                    from './modules.nf'
include { CUSTOM_VARIANT_FILTER_PROVENANCE }                                       from './modules.nf'
include { LIST_SAMPLES_FROM_GROUP_VCF }                                            from './modules.nf'
include { ANNOTATE_SAMPLE_VCF }                                                    from './modules.nf'
include { IDENTIFY_GERMLINE_FROM_STATS }                                           from './modules.nf'
include { EXTRACT_GERMLINE_PREVALENCE_TABLE }                                      from './modules.nf'
include { PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS }                                 from './modules.nf'
include { SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS }                  from './modules.nf'
include { CREATE_ADO_TABLE_FROM_GERMLINE_VCF }                                     from './modules.nf'
include { SUMMARIZE_ADO_INTERVALS }                                                from './modules.nf'
include { PLOT_ADO_GERMLINE_COMPARISON }                                           from './modules.nf'
include { CONCAT_SUMMARY_ADO_INTERVALS_LABELED as CONCAT_ADO_STATS }               from './modules.nf'
include { CONCAT_SUMMARY_ADO_INTERVALS_LABELED as CONCAT_ADO_VEP }                 from './modules.nf'
include { CONCAT_SUMMARY_ADO_INTERVALS_LABELED as CONCAT_ADO_BULK }                from './modules.nf'
// ── New-DAG modules (per-cell ML filter, per-chr VEP, VAF-split cascade, focal
//    pileup, scheme matrices, and QC/filter diagnostics) ──
include { FILTER_MERGED_VCF_BY_ML_VERDICT }                                        from './modules.nf'
include { FILTER_VEP_GERMLINE_CHR }                                                from './modules.nf'
include { FILTER_CHOSEN_VARIANTS_BY_VEP_CHR }                                      from './modules.nf'
include { CONCAT_VEP_PERCHR_OUTPUTS }                                              from './modules.nf'
include { MERGE_MANDATORY_PRIORITY_VARIANTS }                                      from './modules.nf'
include { SEQUOIA_SECOND_FILTER_MERGE }                                            from './modules.nf'
include { VARIANT_PREFILTER_TABLE }                                                from './modules.nf'
include { VARIANT_PREFILTER_TABLE_PERCHR }                                         from './modules.nf'
include { VAF_SPLIT_VARIANTS_HEXBIN }                                              from './modules.nf'
include { VARIANT_FILTER_FUNNEL }                                                  from './modules.nf'
include { BUILD_FOCAL_PILEUP }                                                     from './modules.nf'
include { CONCAT_FOCAL_PILEUP }                                                    from './modules.nf'
include { EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF }                                    from './modules.nf'
include { CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF }                              from './modules.nf'
include { ANNOTATE_VCF_SCHEME_MEMBERSHIP }                                         from './modules.nf'
include { MERGE_ANNOTATED_SAMPLE_VCFS }                                            from './modules.nf'
include { SUBSET_ANNOTATED_VCF_HQSTAT_QC_DEPTH }                                   from './modules.nf'
include { COMPARE_HQSTAT_QC_REDUNDANCY }                                           from './modules.nf'
include { QUANTIFY_QC_FILTER_INFLUENCE }                                           from './modules.nf'
include { MANDATORY_VARIANTS_QC_STATUS }                                           from './modules.nf'
include { COLLECT_DOWNSTREAM_ARTIFACTS }                                           from './modules.nf'
include { EMIT_LINEAGE_INPUTS_MANIFEST }                                           from './modules.nf'
include { EXPLORE_QC_FILTER_THRESHOLDS }                                           from './modules.nf'
include { COHORT_METRICS_SUMMARY }                                                 from './modules.nf'
include { ANALYZE_NR_NV_PILEUP_DEPTH }                                             from './modules.nf'
include { PLOT_MATRIX_SCHEME_SUMMARY }                                             from './modules.nf'

// ============================================================================
// PROCESS: SENTIEON_DNASCOPE  (inline — architecture-specific container)
// Description: Run Sentieon DNAscope germline variant calling on BAM/CRAM → VCF.
//              Two-step: DNAscope (raw calls) + DNAModelApply (ML filter).
//              Supports both BAM and CRAM input formats.
// ============================================================================
process SENTIEON_DNASCOPE {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(alignment_file), path(index_file)
    path(fasta_ref)
    val(interval)
    path(model)
    val(pcrfree)
    val(ploidy)
    val(emit_mode)

    output:
    tuple val(sample_name), path("${sample_name}_dnascope.vcf.gz"), emit: vcf
    path("sentieon_dnascope_${sample_name}_version.yml"),            emit: version

    script:
    def interval_param = interval ? "--interval ${interval}" : ''
    def pcr_indel_model = pcrfree.toString().toBoolean() ? "--pcr_indel_model none" : ''
    """
    set +u
    export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

    sentieon driver \\
        -t ${task.cpus} \\
        -r ${fasta_ref}/genome.fa \\
        -i ${alignment_file} \\
        ${interval_param} \\
        --algo DNAscope \\
        ${pcr_indel_model} \\
        --ploidy ${ploidy} \\
        --model ${model}/dnascope.model \\
        --emit_mode ${emit_mode} \\
        TMP_VARIANT.vcf.gz

    sentieon driver \\
        -t ${task.cpus} \\
        -r ${fasta_ref}/genome.fa \\
        --algo DNAModelApply \\
        --model ${model}/dnascope.model \\
        -v TMP_VARIANT.vcf.gz \\
        ${sample_name}_dnascope.vcf.gz

    export SENTIEON_VER="202503.02"
    echo "Sentieon: \$SENTIEON_VER" > sentieon_dnascope_${sample_name}_version.yml
    """
}

// ============================================================================
// PROCESS: VCF_TO_PARQUET_SOMATIC_VARIANTS  (inline — pipeline-specific)
// Description: Convert per-sample somatic VCF to Parquet with Hive-style
//              partitioning for the BaseJumper lakehouse somatic_variants_summary table.
// Output path: somatic_variants_summary/workspace=*/workflow_id=*/biosample=*/output.parquet
// ============================================================================
process VCF_TO_PARQUET_SOMATIC_VARIANTS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(vcf_files)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)

    output:
    path("somatic_variants_summary/workspace=*/workflow_id=*/biosample=*/output.parquet"), emit: parquet

    script:
    """
    python3 << 'PYEOF'
import gzip, os, glob
from collections import defaultdict
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

workspace        = "${workspace}"
workflow_id      = "${workflow_id}"
biosample        = "${sample_name}"
pipeline_version = "${pipeline_version}"

# Locate the staged vcf.gz (exclude .tbi)
vcf_files = [f for f in glob.glob("*_somatic_annotated.vcf.gz") if not f.endswith(".tbi")]
assert vcf_files, "No *_somatic_annotated.vcf.gz found"
vcf_path = vcf_files[0]

# Columnar approach: append to per-column lists, build DataFrame once at the end.
# Much faster and lower-memory than list-of-dicts for 50k+ variants.
# Handles variable FORMAT fields across rows by backfilling "." for new/missing keys.
data         = defaultdict(list)
sample_names = []
n_rows       = 0

FIXED_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

opener = gzip.open if vcf_path.endswith(".gz") else open

with opener(vcf_path, "rt") as fh:
    for line in fh:
        line = line.rstrip("\\n")
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [sample ...]
            cols         = line.lstrip("#").split("\\t")
            sample_names = cols[9:]
            continue
        # Data row — split once, access by index (faster than dict per row)
        parts = line.split("\\t")
        data["CHROM"].append(parts[0])
        data["POS"].append(int(parts[1]))
        data["ID"].append(parts[2])
        data["REF"].append(parts[3])
        data["ALT"].append(parts[4])
        data["QUAL"].append(parts[5])
        data["FILTER"].append(parts[6])
        fmt = parts[8].split(":")
        cols_set = set(FIXED_COLS)
        for s_idx, s_name in enumerate(sample_names):
            s_vals = parts[9 + s_idx].split(":")
            for k_idx, key in enumerate(fmt):
                col_name = f"{s_name}_{key}" if len(sample_names) > 1 else key
                # First time this FORMAT key is seen: backfill all previous rows with "."
                if col_name not in data:
                    data[col_name] = ["."] * n_rows
                data[col_name].append(s_vals[k_idx] if k_idx < len(s_vals) else ".")
                cols_set.add(col_name)
        # Fill columns that existed before but are absent from this row's FORMAT
        for col_name in list(data.keys()):
            if col_name not in cols_set:
                data[col_name].append(".")
        n_rows += 1

# Single DataFrame construction from column arrays — O(n), no per-row dict overhead
df = pd.DataFrame(data)

# Prepend metadata columns
for i, (col, val) in enumerate([
    ("pipeline_version", pipeline_version),
    ("workspace",        workspace),
    ("workflow_id",      workflow_id),
    ("biosample",        biosample),
]):
    df.insert(i, col, val)

# ── Schema enforcement ──────────────────────────────────────────────────────
_str_cols    = ["pipeline_version", "workspace", "workflow_id", "biosample",
                "CHROM", "ID", "REF", "ALT", "FILTER",
                "GT", "AD", "ML", "PL", "SB", "PGT", "PID"]
_float_cols  = ["QUAL"]
_int_cols    = ["POS", "DP", "GQ"]

for col in _float_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("float64")
for col in _int_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
for col in _str_cols:
    if col in df.columns:
        df[col] = df[col].astype("string")
# ────────────────────────────────────────────────────────────────────────────

out_dir = f"somatic_variants_summary/workspace={workspace}/workflow_id={workflow_id}/biosample={biosample}"
os.makedirs(out_dir, exist_ok=True)

pq.write_table(
    pa.Table.from_pandas(df, preserve_index=False),
    os.path.join(out_dir, "output.parquet")
)
print(f"Written {len(df)} variants to {out_dir}/output.parquet")
PYEOF
    """
}

// ============================================================================
// PROCESS: SOMATIC_STATS_FROM_ANNOTATED_VCF  (inline — pipeline-specific)
// Description: Derive per-sample somatic variant counts from per-sample
//              annotated VCFs produced by ANNOTATE_SAMPLE_VCF.
//              Variants where GT == 0/0 are reference (sample does not carry
//              the variant); all others are counted as somatic.
//              Emits one MultiQC custom-content TSV + group HTML summary.
// ============================================================================
process SOMATIC_STATS_FROM_ANNOTATED_VCF {
    tag "${group}"

    input:
    tuple val(group), path(vcf_files)

    output:
    path("somatic_variant_summary_mqc.tsv"), emit: stats
    path("somatic_group_summary_mqc.html"),  emit: summary

    script:
    """
    python3 << 'PYEOF'
import gzip, glob, statistics

vcf_paths    = sorted([f for f in glob.glob("*_somatic_annotated.vcf.gz")
                       if not f.endswith(".tbi")])
sample_names = [f.replace("_somatic_annotated.vcf.gz", "") for f in vcf_paths]

counts            = {}
n_cohort_variants = 0

for vcf_path, sample in zip(vcf_paths, sample_names):
    n_total = n_snv = n_indel = 0
    n_variants_in_vcf = 0
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\\n").split("\\t")
            ref, alt = parts[3], parts[4]
            n_variants_in_vcf += 1

            # Parse GT to determine if sample carries this variant.
            # ANNOTATE_SAMPLE_VCF sets missing positions to 0/0 via +setGT.
            fmt_keys = parts[8].split(":")
            smp_vals = parts[9].split(":")
            gt_idx   = fmt_keys.index("GT") if "GT" in fmt_keys else 0
            gt       = smp_vals[gt_idx] if gt_idx < len(smp_vals) else "."
            gt_norm  = gt.replace("|", "/")

            if gt_norm not in ("0/0", ".", "./.", "0"):
                n_total += 1
                if len(ref) != len(alt):
                    n_indel += 1
                else:
                    n_snv += 1

    counts[sample]    = {"n_total": n_total, "n_snv": n_snv, "n_indel": n_indel}
    n_cohort_variants = max(n_cohort_variants, n_variants_in_vcf)

with open("somatic_variant_summary_mqc.tsv", "w") as out:
    out.write("# id: 'somatic_variant_summary'\\n")
    out.write("# section_name: 'Somatic Summary'\\n")
    out.write("# plot_type: 'table'\\n")
    out.write("# headers:\\n")
    out.write("#   n_total:\\n")
    out.write("#     title: 'Total Somatic'\\n")
    out.write("#     scale: 'Blues'\\n")
    out.write("#   n_snv:\\n")
    out.write("#     title: 'SNVs'\\n")
    out.write("#     scale: 'Greens'\\n")
    out.write("#   n_indel:\\n")
    out.write("#     title: 'INDELs'\\n")
    out.write("#     scale: 'Oranges'\\n")
    out.write("Sample\\tn_total\\tn_snv\\tn_indel\\n")
    for sample in sample_names:
        c = counts[sample]
        out.write(f"{sample}\\t{c['n_total']}\\t{c['n_snv']}\\t{c['n_indel']}\\n")

totals     = [counts[s]["n_total"] for s in sample_names]
n_cells    = len(totals)
n_with_var = sum(1 for t in totals if t > 0)
med        = statistics.median(totals) if totals else 0
avg        = sum(totals) / n_cells if n_cells else 0

def row(label, val):
    return (f"<tr><td>{label}</td>"
            f"<td style='text-align:right'>{val}</td></tr>\\n")

with open("somatic_group_summary_mqc.html", "w") as out:
    out.write("<!--\\nid: 'somatic_group_summary'\\nsection_name: 'Somatic Cohort Summary'\\nplot_type: 'html'\\n-->\\n")
    out.write("<div style='padding:12px'>\\n")
    out.write("<table class='table table-bordered table-condensed mqc-table' style='min-width:360px'>\\n")
    out.write("<thead><tr><th>Metric</th><th style='text-align:right'>Value</th></tr></thead>\\n")
    out.write("<tbody>\\n")
    out.write(row("Total cohort variant sites (annotated VCF)", f"{n_cohort_variants:,}"))
    out.write(row("Cells analyzed", f"{n_cells:,}"))
    out.write(row("Cells with &ge;1 somatic variant", f"{n_with_var:,}"))
    out.write(row("Min variants per cell",    f"{min(totals):,}" if totals else "0"))
    out.write(row("Median variants per cell", f"{med:.1f}"))
    out.write(row("Mean variants per cell",   f"{avg:.1f}"))
    out.write(row("Max variants per cell",    f"{max(totals):,}" if totals else "0"))
    out.write("</tbody>\\n</table>\\n</div>\\n")
PYEOF
    """
}

// ============================================================================
// PROCESS: MULTIQC_SOMATIC  (inline — pipeline-specific)
// Description: Aggregate per-sample somatic stats TSVs and (when VEP is
//              enabled) per-chromosome VEP summary HTMLs into a single
//              MultiQC HTML report.
// ============================================================================
process MULTIQC_SOMATIC {
    input:
    path(input_files)
    val(project)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    path(logo)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    """
    cat > multiqc_config.yaml << EOF
custom_logo_title: 'BioSkryb Genomics'
custom_logo: bioskryb_logo-tagline.png
custom_logo_width: 260

title: "basej-somatic v${pipeline_version}"
report_header_info:
  - Project: "${project}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false

run_modules:
  - custom_content
  - vep

custom_content:
  order:
    - somatic_variant_summary
    - somatic_group_summary
    - somatic_ado_analysis
    - somatic_filter_funnel
    - somatic_matrix_scheme
    - somatic_hexbin_first
    - somatic_hexbin_second
    - somatic_threshold_shared
    - somatic_threshold_singleton
    - somatic_hqstat_qc_redundancy
    - somatic_pileup_sparsity
    - somatic_qc_influence
    - somatic_qc_threshold_as_propclipped
    - somatic_qc_threshold_bppos
    - somatic_qc_threshold_passrates
    - somatic_germline_prevalence
    - somatic_germline_plots
    - somatic_digitalheatmap
    - somatic_vafheatmap
    - somatic_provenance_report

EOF

    python3 << 'PYEOF'
import base64, glob

# ── June 2026: filtering-diagnostics panels (parity with Isai's MULTIQC_REPORT) ──
# Helper: embed one or more PNGs as an HTML custom-content section (our format).
def _embed_pngs(section_id, title, png_glob, description=""):
    pngs = sorted(glob.glob(png_glob))
    if not pngs:
        return
    with open(f"{section_id}_mqc.html", "w") as out:
        out.write(f"<!--\\nid: '{section_id}'\\nsection_name: '{title}'\\nplot_type: 'html'\\n-->\\n")
        out.write("<div style='padding:12px'>\\n")
        if description:
            out.write(f"<p>{description}</p>\\n")
        for p in pngs:
            with open(p, "rb") as f:
                b64 = base64.b64encode(f.read()).decode()
            out.write(f"<img src='data:image/png;base64,{b64}' style='width:100%; margin-bottom:12px; display:block;' />\\n")
        out.write("</div>\\n")

# Cohort filter metrics section intentionally removed (not needed in the report).

_embed_pngs("somatic_filter_funnel", "Variant Filter Funnel",
            "variant_filter_plot_*.png",
            "Per-stage variant counts through the full somatic filter cascade.")
_embed_pngs("somatic_matrix_scheme", "Matrix Scheme Summary",
            "matrix_scheme_summary_*_page-*.png",
            "SNV/indel counts per NR/NV filtering scheme (cohort and per-cell).")
_embed_pngs("somatic_hexbin_first", "Rho vs Germline q-value — First-round Sequoia",
            "*_hexbin_FirstRunSequoia_Rho_vs_GermlineQval_all.png")
_embed_pngs("somatic_hexbin_second", "Rho vs Germline q-value — Second-round Sequoia (shared)",
            "*_hexbin_SecondRunSequoia_Rho_vs_GermlineQval_shared.png")
_embed_pngs("somatic_threshold_shared", "VAF-split Threshold Exploration — shared variants",
            "*_threshold_exploration_filtered_shared.png")
_embed_pngs("somatic_threshold_singleton", "VAF-split Threshold Exploration — singleton variants",
            "*_threshold_exploration_filtered_singleton.png")
_embed_pngs("somatic_hqstat_qc_redundancy", "HQStat vs Pileup-QC Redundancy",
            "hqstat_vs_qc_heatmap_*.png",
            "Is the HQ-statistical filter redundant with the pileup-QC filter?")
_embed_pngs("somatic_pileup_sparsity", "Pileup VAF Sparsity",
            "pileup_sparsity_hexbin_*_page-*.png")
_embed_pngs("somatic_qc_influence", "QC Filter Influence",
            "qc_filter_influence_*.png",
            "Influence of each artifact-QC filter (AS / PropClipped / BPPos) on variant removal.")
_embed_pngs("somatic_qc_threshold_as_propclipped", "QC Threshold Calibration — AS / PropClipped",
            "qc_threshold_as_propclipped_*.png")
_embed_pngs("somatic_qc_threshold_bppos", "QC Threshold Calibration — BP-position",
            "qc_threshold_bppos_*.png")
_embed_pngs("somatic_qc_threshold_passrates", "QC Threshold Calibration — Pass-rate sensitivity",
            "qc_threshold_passrates_*.png")

# Somatic variant summary & group summary tables (from SOMATIC_STATS_FROM_ANNOTATED_VCF)
# These will be picked up by MultiQC if they follow the naming convention

# Small single-page PDFs (digitalheatmap, vafheatmap) embedded directly
pdf_meta = {
    "digitalheatmap": ("somatic_digitalheatmap", "Digital Heatmap"),
    "vafheatmap":     ("somatic_vafheatmap",      "VAF Heatmap"),
}
for pdf_path in sorted(glob.glob("*.pdf")):
    for key, (section_id, title) in pdf_meta.items():
        if key in pdf_path:
            with open(pdf_path, "rb") as f:
                b64 = base64.b64encode(f.read()).decode()
            with open(f"{section_id}_mqc.html", "w") as out:
                out.write(f"<!--\\nid: '{section_id}'\\nsection_name: '{title}'\\nplot_type: 'html'\\n-->\\n")
                out.write("<div style='padding:12px'>\\n")
                out.write(f"<embed type='application/pdf' src='data:application/pdf;base64,{b64}' width='100%' height='800px' />\\n")
                out.write("</div>\\n")
            break

# Germline prevalence table — TSV with per-variant filter statistics
prevalence_tsvs = glob.glob("germline_prevalence_long_*.tsv")
if prevalence_tsvs:
    with open(prevalence_tsvs[0]) as f:
        content = f.read()
    with open("somatic_germline_prevalence_mqc.tsv", "w") as out:
        out.write("# id: 'somatic_germline_prevalence'\\n")
        out.write("# section_name: 'Germline Prevalence Table'\\n")
        out.write("# plot_type: 'table'\\n")
        out.write("# description: 'Per-variant germline prevalence across three filter sets (GERMLINE_FROM_STATS, VEP_AF_filter, Bulk_Fail)'\\n")
        out.write("# pconfig:\\n")
        out.write("#   id: 'germline_prevalence_table'\\n")
        out.write(content)

# Germline prevalence distribution plots — PNG embedded as image
prevalence_pngs = glob.glob("germline_prevalence_distributions_*.png")
if prevalence_pngs:
    with open("somatic_germline_plots_mqc.html", "w") as out:
        out.write("<!--\\nid: 'somatic_germline_plots'\\nsection_name: 'Germline Prevalence Distributions'\\nplot_type: 'html'\\n-->\\n")
        out.write("<div style='padding:12px'>\\n")
        out.write("<p>Distribution of germline variant prevalence across cell populations, stratified by filter set.</p>\\n")
        for png_path in sorted(prevalence_pngs):
            with open(png_path, "rb") as f:
                b64 = base64.b64encode(f.read()).decode()
            group_name = png_path.replace("germline_prevalence_distributions_", "").replace(".png", "")
            out.write(f"<h4>Group: {group_name}</h4>\\n")
            out.write(f"<img src='data:image/png;base64,{b64}' style='width:100%; margin-bottom:16px; display:block;' />\\n")
        out.write("</div>\\n")

# ADO (Allele Drop-Out) analysis results — PNG plots embedded as images
# ADO (Allele Drop-Out) analysis results — PNG plots embedded as images.
# The "Allelic balance: Summary" plot (ADO_germline_summary.png) is intentionally
# excluded from the report; comparison + distribution plots are kept.
ado_plots = [p for p in sorted(glob.glob("ADO_germline_*.png"))
             if not p.endswith("ADO_germline_summary.png")]
ado_summary_tsv = glob.glob("ADO_germline_summary.tsv")

if ado_plots or ado_summary_tsv:
    with open("somatic_ado_analysis_mqc.html", "w") as out:
        out.write("<!--\\nid: 'somatic_ado_analysis'\\nsection_name: 'ADO (Allele Drop-Out) Analysis'\\nplot_type: 'html'\\n-->\\n")
        out.write("<div style='padding:12px'>\\n")
        out.write("<p>Allele drop-out (ADO) rates and patterns across samples, compared across three germline filter sets (stats/vep/bulk).</p>\\n")
        
        # Embed ADO comparison plots
        for png_path in ado_plots:
            with open(png_path, "rb") as f:
                b64 = base64.b64encode(f.read()).decode()
            plot_name = png_path.replace("ADO_germline_", "").replace(".png", "").replace("_", " ").title()
            out.write(f"<h4>{plot_name}</h4>\\n")
            out.write(f"<img src='data:image/png;base64,{b64}' style='width:100%; margin-bottom:16px; display:block;' />\\n")
        
        # Embed ADO summary table if available
        if ado_summary_tsv:
            with open(ado_summary_tsv[0]) as f:
                lines = f.readlines()
            out.write("<h4>ADO Summary Statistics</h4>\\n")
            out.write("<table style='font-family:monospace; font-size:11px; border-collapse:collapse;'>\\n")
            for i, line in enumerate(lines[:20]):  # Limit to first 20 rows for display
                if i == 0:
                    out.write("<thead><tr>")
                    for col in line.rstrip().split("\\t"):
                        out.write(f"<th style='border:1px solid #ccc; padding:4px;'>{col}</th>")
                    out.write("</tr></thead><tbody>\\n")
                else:
                    out.write("<tr>")
                    for col in line.rstrip().split("\\t"):
                        out.write(f"<td style='border:1px solid #ccc; padding:4px;'>{col}</td>")
                    out.write("</tr>\\n")
            out.write("</tbody></table>\\n")
        
        out.write("</div>\\n")

PYEOF

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

    if (!params.platform) {
        exit 1, "ERROR: --platform parameter is required (Illumina or Ultima)"
    }

    // Chr list is derived from gender and genome:
    //   GRCm39: chr1-chr19 (19 autosomes) + chrX/Y (males get Y, females only X)
    //   GRCh38: chr1-chr22 (22 autosomes) + chrX/Y (males get Y, females only X)
    // Can be overridden with --chrs on the CLI for targeted runs.
    def default_chrs_male_grch38   = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    def default_chrs_female_grch38 = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
    def default_chrs_male_grcm39   = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
    def default_chrs_female_grcm39 = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX']
    
    def default_chrs_male = params.genome == "GRCm39" ? default_chrs_male_grcm39 : default_chrs_male_grch38
    def default_chrs_female = params.genome == "GRCm39" ? default_chrs_female_grcm39 : default_chrs_female_grch38
    def chrs_list = params.chrs
        ? ((params.chrs instanceof List) ? params.chrs : params.chrs.replaceAll(/[\[\]\s"']/, '').split(',').toList())
        : (params.gender == 'male' ? default_chrs_male : default_chrs_female)

    // -------------------------------------------------------------------------
    // Phase 0: Parse input CSV (platform-aware)
    // 
    // Illumina:
    //   Required columns: biosampleName, bam, group
    //   Index (bai) auto-discovered from bam + ".bai"
    //   Optional column: vcf — pre-computed DNAscope VCF (.vcf.gz) to skip calling
    //
    // Ultima:
    //   Required columns: biosampleName, cram, group
    //   Index (crai) auto-discovered from cram + ".crai"
    //   Optional column: vcf — pre-computed DNAscope VCF (.vcf.gz) to skip calling
    // -------------------------------------------------------------------------
    ch_input = Channel.fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (params.platform.equalsIgnoreCase("Ultima")) {
                // Ultima: CRAM input
                def cram = file(row.cram)
                def crai = row.crai ? file(row.crai) : file(row.cram + ".crai")
                tuple(
                    row.biosampleName,
                    cram,
                    crai,
                    row.groups ?: row.group ?: 'GROUP1',
                    row.vcf ?: ''
                )
            } else if (params.platform.equalsIgnoreCase("Illumina") || params.platform.equalsIgnoreCase("Element")) {
                // Illumina/Element: BAM input
                tuple(
                    row.biosampleName,
                    file(row.bam),
                    file(row.bam + ".bai"),
                    row.groups ?: row.group ?: 'GROUP1',
                    row.vcf ?: ''
                )
            } else {
                exit 1, "ERROR: --platform must be either 'Illumina', 'Element', or 'Ultima'"
            }
        }

    ch_input.ifEmpty { exit 1, "ERROR: No samples found in --input_csv" }

    // -------------------------------------------------------------------------
    // Phase 1: BAM → VCF via Sentieon DNAscope (skipped when vcf column is set)
    // -------------------------------------------------------------------------

    // Branch on the vcf field of the channel tuple (not CSV column position;
    // CSV columns are accessed by name so their order does not matter).
    // Both BAM and CRAM flow through the same logic.
    ch_input.branch { name, alignment_file, index_file, group, vcf ->
        has_vcf:        vcf != ''
        needs_dnascope: true
    }.set { ch_branched }

    SENTIEON_DNASCOPE(
        ch_branched.needs_dnascope.map { name, alignment_file, index_file, group, vcf -> tuple(name, alignment_file, index_file) },
        params.reference,
        params.calling_intervals,
        params.dnascope_model,
        params.pcrfree,
        params.ploidy,
        params.dnascope_emit_mode
    )

    // -------------------------------------------------------------------------
    // Build input channels for somatic filtering
    // ch_alignment: (biosampleName, [alignment_file, index_file], group) — for pileup
    //               (works with both BAM and CRAM)
    // ch_vcf: (biosampleName, [vcf_gz], group)      — for PREPROCESS_VCF
    // -------------------------------------------------------------------------
    ch_alignment = ch_input.map { name, alignment_file, index_file, group, vcf -> tuple(name, tuple(alignment_file, index_file), group) }

    ch_group_map = ch_input.map { name, alignment_file, index_file, group, vcf -> tuple(name, group) }

    ch_vcf_from_dnascope = SENTIEON_DNASCOPE.out.vcf
        .map { name, vcf -> tuple(name, [vcf]) }
        .join(ch_group_map)
        .map { name, vcf_list, group -> tuple(name, vcf_list, group) }

    ch_vcf_precomputed = ch_branched.has_vcf
        .map { name, alignment_file, index_file, group, vcf -> tuple(name, [file(vcf)]) }
        .join(ch_group_map)
        .map { name, vcf_list, group -> tuple(name, vcf_list, group) }

    ch_vcf = ch_vcf_from_dnascope.mix(ch_vcf_precomputed)

    // -------------------------------------------------------------------------
    // Phase 2: Preprocess and merge per-sample VCFs
    // -------------------------------------------------------------------------
    PREPROCESS_VCF(
        ch_vcf,
        params.reference,
        params.model_vcf
    )

    ch_input_merge_processed_vcf = PREPROCESS_VCF.out.vcf
        .map { sample_name, vcf_files, group -> tuple(group, vcf_files) }
        .groupTuple(by: 0)
        .map { group, vcf_lists -> tuple(group, vcf_lists.flatten().collect()) }

    MERGE_PROCESSED_VCF(
        ch_input_merge_processed_vcf,
        params.reference
    )

    // -------------------------------------------------------------------------
    // Phase 3: Per-cell ML-verdict filter — subset merged VCF to "usable" variants
    // MERGE_PROCESSED_VCF carries each cell's caller verdict as FORMAT/MLV. This
    // module keeps variants accepted by >= ml_min_pass_cells cells. All downstream
    // consumes the usable VCF + its regenerated df_nv.
    // -------------------------------------------------------------------------
    FILTER_MERGED_VCF_BY_ML_VERDICT(
        MERGE_PROCESSED_VCF.out.merged_vcf,
        params.ml_min_pass_cells
    )

    ch_usable_merged_vcf = FILTER_MERGED_VCF_BY_ML_VERDICT.out.usable_vcf
    ch_usable_df_nv      = FILTER_MERGED_VCF_BY_ML_VERDICT.out.df_nv

    // Unpacked (group, vcf, tbi) form of the usable merged VCF for modules that take
    // vcf + tbi as separate inputs (the .out emits a glob list [vcf, tbi]).
    ch_usable_vcf_unpacked = ch_usable_merged_vcf
        .map { group, files ->
            def vcf = files.find { it.name.endsWith('.vcf.gz') && !it.name.endsWith('.tbi') }
            def tbi = files.find { it.name.endsWith('.tbi') }
            tuple(group, vcf, tbi)
        }

    // -------------------------------------------------------------------------
    // Phase 4: Bulk variant list (optional)
    // -------------------------------------------------------------------------
    if (params.bulk_vcf != "") {
        ch_bulk_vcf = Channel.fromPath(params.bulk_vcf)
        BULK_GET_VARIANTS_TO_FILTER(
            ch_bulk_vcf,
            params.reference,
            params.model_vcf
        )
        ch_bulk_variants_to_remove = BULK_GET_VARIANTS_TO_FILTER.out
    } else {
        CREATE_EMPTY_BULK_VARIANTS(Channel.value(true))
        ch_bulk_variants_to_remove = CREATE_EMPTY_BULK_VARIANTS.out.bulk_variants
    }

    GET_VARIANTS_FROM_MERGED_VCF(
        ch_usable_merged_vcf
    )

    // Effective chromosome list (gender/genome-aware chrs_list + optional extras).
    def chrs_eff = (params.extra_chromosomes && params.extra_chromosomes != "")
        ? chrs_list + params.extra_chromosomes.tokenize(',').collect { it.trim() }
        : chrs_list

    // -------------------------------------------------------------------------
    // Phase 5: VEP per-chromosome germline filtering (runs BEFORE the bulk filter)
    // When vep_cache_dir is empty the VEP processes are skipped and the chosen
    // list falls back to the full GET_VARIANTS_FROM_MERGED_VCF output.
    // -------------------------------------------------------------------------
    if (params.vep_cache_dir) {
        ch_vep_cache = file(params.vep_cache_dir, type: 'dir')

        ch_merged_by_chr = ch_usable_merged_vcf
            .combine(Channel.of(chrs_eff).flatMap())
            .map { group, vcf_tbi_list, chr -> tuple(group, vcf_tbi_list[0], vcf_tbi_list[1], chr) }

        SPLIT_SUBSET_VCF_BY_CHR(
            ch_merged_by_chr
        )

        VEP_ANNOTATE(
            SPLIT_SUBSET_VCF_BY_CHR.out.subset_vcf_chr,
            params.reference,
            params.vep_species,
            params.vep_assembly,
            ch_vep_cache
        )

        SORT_INDEX_VEP_VCF(
            VEP_ANNOTATE.out.vep_vcf_chr
        )

        // Group-level VEP VCF (for per-sample annotation) — sorted by chr order.
        ch_vep_by_group = SORT_INDEX_VEP_VCF.out.vep_vcf_chr
            .groupTuple(by: 0)
            .map { group, list_chr, list_vcf, list_tbi ->
                tuple(group, [list_chr, list_vcf].transpose().sort { a, b -> a[0] <=> b[0] }.collect { it[1] })
            }

        MERGE_VEP_VCF_BY_GROUP(
            ch_vep_by_group
        )

        // Per-chromosome germline filter (parallel; no wait for the merge).
        ch_filter_vep_germline_chr_in = SORT_INDEX_VEP_VCF.out.vep_vcf_chr
            .map { group, chr, vcf, tbi -> tuple(group, chr, vcf) }

        FILTER_VEP_GERMLINE_CHR(
            ch_filter_vep_germline_chr_in,
            params.vep_max_af_1kg,
            params.vep_filter_by_existing_variation
        )

        // Intersect each chr's VEP-germline-pass list with the full all-variants list.
        ch_filter_vep_input = FILTER_VEP_GERMLINE_CHR.out.chosen_variants
            .combine(GET_VARIANTS_FROM_MERGED_VCF.out.all_variants, by: 0)
            .map { group, chr, vep_pass_chr, all_variants -> tuple(group, chr, all_variants, vep_pass_chr) }

        FILTER_CHOSEN_VARIANTS_BY_VEP_CHR(
            ch_filter_vep_input
        )

        // Reassemble per-chr outputs into group-level files (one CONCAT per group).
        ch_chosen_chr_by_group = FILTER_CHOSEN_VARIANTS_BY_VEP_CHR.out.chosen_variants
            .map { group, chr, f -> tuple(group, f) }
            .groupTuple(by: 0)
        ch_provenance_chr_by_group = FILTER_VEP_GERMLINE_CHR.out.filter_provenance
            .map { group, chr, f -> tuple(group, f) }
            .groupTuple(by: 0)
        ch_priority_chr_by_group = FILTER_VEP_GERMLINE_CHR.out.priority_variants
            .map { group, chr, f -> tuple(group, f) }
            .groupTuple(by: 0)

        ch_concat_vep_input = ch_chosen_chr_by_group
            .join(ch_provenance_chr_by_group, by: 0)
            .join(ch_priority_chr_by_group, by: 0)

        CONCAT_VEP_PERCHR_OUTPUTS(
            ch_concat_vep_input
        )

        ch_chosen_pre_bulk            = CONCAT_VEP_PERCHR_OUTPUTS.out.chosen_variants
        ch_vep_filter_provenance      = CONCAT_VEP_PERCHR_OUTPUTS.out.filter_provenance
        ch_merged_vep_vcf_for_samples = MERGE_VEP_VCF_BY_GROUP.out.vep_vcf
            .map { group, files ->
                def vcf = files.find { it.name.endsWith('.vcf.gz') && !it.name.endsWith('.tbi') }
                def tbi = files.find { it.name.endsWith('.vcf.gz.tbi') }
                tuple(group, vcf, tbi)
            }
        ch_priority_variants          = CONCAT_VEP_PERCHR_OUTPUTS.out.priority_variants
        // VEP per-chromosome summary HTMLs (for MultiQC) + merged VEP VCF (for publish)
        ch_vep_html_for_multiqc       = VEP_ANNOTATE.out.vep_stats.map { group, chr, html -> html }
        ch_vep_vcf_for_publish        = MERGE_VEP_VCF_BY_GROUP.out.vep_vcf
    } else {
        // VEP disabled: pass the full all-variants list straight to the bulk filter,
        // and substitute the usable merged VCF + sentinels for the VEP-derived channels.
        ch_chosen_pre_bulk            = GET_VARIANTS_FROM_MERGED_VCF.out.all_variants
        ch_vep_filter_provenance      = ch_usable_merged_vcf.map { group, _l -> tuple(group, file("${projectDir}/assets/dummy_file.txt")) }
        ch_merged_vep_vcf_for_samples = ch_usable_vcf_unpacked
        ch_priority_variants          = ch_usable_merged_vcf.map { group, _l -> tuple(group, file("${projectDir}/assets/dummy_file.txt")) }
        // No VEP: empty channels so MultiQC and the vep_annotated publish are skipped.
        ch_vep_html_for_multiqc       = Channel.empty()
        ch_vep_vcf_for_publish        = Channel.empty()
    }

    // -------------------------------------------------------------------------
    // Phase 6: Merge mandatory variants into the priority list (bypass all filters)
    // -------------------------------------------------------------------------
    if (params.mandatory_variants) {
        ch_mandatory_variants_file = Channel.fromPath(params.mandatory_variants)
    } else {
        ch_mandatory_variants_file = Channel.of(file("${projectDir}/assets/dummy_file.txt"))
    }

    MERGE_MANDATORY_PRIORITY_VARIANTS(
        ch_priority_variants.combine(ch_mandatory_variants_file)
    )

    ch_priority_variants = MERGE_MANDATORY_PRIORITY_VARIANTS.out.merged_priority_variants

    // -------------------------------------------------------------------------
    // Phase 7: Bulk filter on the (VEP-filtered) chosen variants
    // Priority/mandatory variants (ch_priority_variants) are rescued: a variant in
    // the bulk list that is also priority is kept and marked "Rescued" in provenance.
    // -------------------------------------------------------------------------
    ch_filter_bulk_input_early = ch_chosen_pre_bulk
        .combine(ch_bulk_variants_to_remove)
        .combine(ch_priority_variants, by: 0)

    FILTER_CHOSEN_VARIANTS_BY_BULK(
        ch_filter_bulk_input_early
    )

    ch_chosen_after_bulk_early = FILTER_CHOSEN_VARIANTS_BY_BULK.out.chosen_variants

    GET_LIST_POS_FROM_CHOSEN_VARIANTS(
        ch_chosen_after_bulk_early
    )

    // -------------------------------------------------------------------------
    // Phase 8: Per-sample pileup at candidate variant positions
    // -------------------------------------------------------------------------
    ch_chr = Channel.of(chrs_eff).flatMap()

    ch_input_bam_group_pileup = ch_alignment
        .map { name, alignment_pair, group -> tuple(group, name, alignment_pair) }
        .combine(GET_LIST_POS_FROM_CHOSEN_VARIANTS.out.list_pos, by: 0)
        .combine(ch_chr)

    CUSTOM_BAM_GROUP_PILEUP(
        ch_input_bam_group_pileup,
        params.reference
    )

    // -------------------------------------------------------------------------
    // Phase 9: Filter df_nv (usable) by chosen variants; build NV/NR tables
    // -------------------------------------------------------------------------
    ch_input_df_nr = CUSTOM_BAM_GROUP_PILEUP.out.df_nr
        .map { group, chr, sample_name, df_nr -> tuple(group, chr, df_nr) }
        .groupTuple(by: [0, 1])
        .map { group, chr, files -> tuple(group, chr, files.flatten().collect()) }

    ch_input_filter_df_nv = ch_usable_df_nv
        .combine(ch_chosen_after_bulk_early, by: 0)

    FILTER_DF_NV_BY_CHOSEN_VARIANTS(
        ch_input_filter_df_nv
    )

    ch_input_create_tab_nvnr = FILTER_DF_NV_BY_CHOSEN_VARIANTS.out.df_nv_filtered
        .combine(ch_input_df_nr, by: 0)
        .map { group, df_nv, chr, df_nr_files ->
            tuple(group, chr, df_nv, df_nr_files.flatten().collect())
        }

    CREATE_TAB_NVNR(
        ch_input_create_tab_nvnr
    )

    // -------------------------------------------------------------------------
    // Phase 10: First-pass binomial / betabinomial filter (Sequoia)
    // -------------------------------------------------------------------------
    SEQUOIA_BINOM_BETABINOM_TAB_NV_NR(
        CREATE_TAB_NVNR.out,
        params.aggregated_min_mean_depth,
        params.aggregated_max_mean_depth,
        params.gender
    )

    ch_input_concat_filter_bb = SEQUOIA_BINOM_BETABINOM_TAB_NV_NR.out.df_filter
        .map { group, chr, df -> tuple(group, df) }
        .groupTuple(by: 0)
        .map { group, files -> tuple(group, files.flatten().collect()) }
        .combine(ch_priority_variants, by: 0)

    CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR(
        ch_input_concat_filter_bb,
        params.first_pass_binomial_cutoff,
        params.first_pass_betabinomial_cutoff
    )

    ch_chosen_variants = CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.chosen_variants

    // -------------------------------------------------------------------------
    // seq_mode preset resolution (explicit params override the preset).
    // bin_nr is the depth gate (wgs depth-limited; wes deep on-target). disable_bppos
    // defaults TRUE for wes (capture read-start clustering breaks the positional-bias
    // assumption) and FALSE for wgs. Other gates are biology-driven and mode-independent.
    // -------------------------------------------------------------------------
    def VAF_SPLIT_PRESETS = [
        wgs: [ rho_thr: 0.2, log10q_thr: -5, bin_nr: 5,  bin_nv: 3, bin_nv_shared: 2,
               bin_vaf: 0.10, bin_vaf_shared_anchor: 0.30, bin_vaf_singleton: 0.30, disable_bppos: 'FALSE' ],
        wes: [ rho_thr: 0.2, log10q_thr: -5, bin_nr: 10, bin_nv: 3, bin_nv_shared: 2,
               bin_vaf: 0.10, bin_vaf_shared_anchor: 0.30, bin_vaf_singleton: 0.30, disable_bppos: 'TRUE'  ]
    ]
    if (!VAF_SPLIT_PRESETS.containsKey(params.seq_mode))
        error "params.seq_mode must be 'wgs' or 'wes' (got '${params.seq_mode}')"
    def _vp = VAF_SPLIT_PRESETS[params.seq_mode]
    def vsRhoThr       = params.vaf_split_rho_thr               != null ? params.vaf_split_rho_thr               : _vp.rho_thr
    def vsLog10qThr    = params.vaf_split_log10q_thr            != null ? params.vaf_split_log10q_thr            : _vp.log10q_thr
    def vsBinNr        = params.vaf_split_bin_nr                != null ? params.vaf_split_bin_nr                : _vp.bin_nr
    def vsBinNv        = params.vaf_split_bin_nv                != null ? params.vaf_split_bin_nv                : _vp.bin_nv
    def vsBinNvShared  = params.vaf_split_bin_nv_shared         != null ? params.vaf_split_bin_nv_shared         : _vp.bin_nv_shared
    def vsBinVaf       = params.vaf_split_bin_vaf               != null ? params.vaf_split_bin_vaf               : _vp.bin_vaf
    def vsBinVafAnchor = params.vaf_split_bin_vaf_shared_anchor != null ? params.vaf_split_bin_vaf_shared_anchor : _vp.bin_vaf_shared_anchor
    def vsBinVafSingle = params.vaf_split_bin_vaf_singleton     != null ? params.vaf_split_bin_vaf_singleton     : _vp.bin_vaf_singleton
    def vsDisableBppos = params.disable_bppos                   != null ? params.disable_bppos                   : _vp.disable_bppos
    log.info "[seq_mode=${params.seq_mode}] gates -> bin_nr=${vsBinNr} rho_thr=${vsRhoThr} log10q_thr=${vsLog10qThr} bin_nv=${vsBinNv} bin_nv_shared=${vsBinNvShared} bin_vaf=${vsBinVaf} bin_vaf_shared_anchor=${vsBinVafAnchor} bin_vaf_singleton=${vsBinVafSingle} disable_bppos=${vsDisableBppos}"

    // -------------------------------------------------------------------------
    // Phase 11: Heuristic sample-level pileup filtering + raw-table QC
    // -------------------------------------------------------------------------
    ch_input_rscript_filter_1 = CUSTOM_BAM_GROUP_PILEUP.out.pileup
        .combine(ch_chosen_variants, by: 0)

    CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR(
        ch_input_rscript_filter_1,
        params.cutoff_mq_hq,
        params.cutoff_bq_hq,
        params.cutoff_bps_start,
        params.num_lines_read_pileup,
        params.read_length
    )

    ch_input_filter_tables = CUSTOM_RSCRIPT_SOMATICSNP_FILTER_1_SAMPLELEVEL_PROCESS_PILEUP_SAMPLE_CIGAR.out
        .groupTuple(by: [0, 1])
        .map { group, chr, files -> tuple(group, chr, files.flatten().collect()) }
        .combine(ch_chosen_variants, by: 0)
        .combine(ch_priority_variants, by: 0)

    CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES(
        ch_input_filter_tables,
        params.cutoff_as,
        params.cutoff_prop_clipped_reads,
        params.cutoff_prop_bp_under,
        params.cutoff_prop_bp_upper,
        params.cutoff_sd_indiv,
        params.cutoff_mad_indiv,
        params.cutoff_sd_both,
        params.cutoff_mad_both,
        params.cutoff_sd_extreme,
        params.cutoff_mad_extreme,
        params.disable_qc,
        vsDisableBppos
    )

    // -------------------------------------------------------------------------
    // Phase 12: Group-level table assembly + per-chr second-pass Sequoia + merge
    // -------------------------------------------------------------------------
    ch_input_tabs_group_level = CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.tabs_all
        .map { group, chr, mat_nv, mat_nr -> tuple(group, mat_nv, mat_nr) }
        .groupTuple(by: 0)
        .map { group, mat_nv_list, mat_nr_list ->
            tuple(group, mat_nv_list.flatten().collect(), mat_nr_list.flatten().collect())
        }

    CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS(
        ch_input_tabs_group_level
    )

    ch_input_sequoia_second_filter = CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS.out.tabs
        .combine(Channel.of(chrs_eff).flatMap())
        .map { group, mat_nv, mat_nr, chr -> tuple(group, chr, mat_nv, mat_nr) }

    SEQUOIA_SECOND_FILTER(
        ch_input_sequoia_second_filter,
        params.reference,
        params.second_pass_binomial_cutoff,
        params.second_pass_betabinomial_cutoff_rho_snp,
        params.second_pass_betabinomial_cutoff_rho_indel,
        params.aggregated_hq_min_mean_depth,
        params.aggregated_hq_max_mean_depth,
        params.gender,
        params.second_pass_beta_binom_shared
    )

    ch_sequoia_second_filter_by_group = SEQUOIA_SECOND_FILTER.out.df_filter
        .map { group, chr, f -> tuple(group, f) }
        .groupTuple(by: 0)
        .map { group, files -> tuple(group, files.flatten()) }

    SEQUOIA_SECOND_FILTER_MERGE(
        ch_sequoia_second_filter_by_group
    )

    // -------------------------------------------------------------------------
    // Phase 13: Prefilter table (per-chr agg + merge) → VAF-split cascade
    // -------------------------------------------------------------------------
    ch_pileup_per_group = CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.pileup
        .map { group, chr, f -> tuple(group, f) }
        .groupTuple()

    VARIANT_PREFILTER_TABLE_PERCHR(
        CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.pileup
    )

    ch_prefilter_merge_input = VARIANT_PREFILTER_TABLE_PERCHR.out.agg
        .groupTuple(by: 0)
        .combine(SEQUOIA_SECOND_FILTER_MERGE.out.df_filter, by: 0)

    VARIANT_PREFILTER_TABLE(
        ch_prefilter_merge_input
    )

    ch_vaf_split_input = CUSTOM_CREATE_GROUP_LEVEL_TAB_DFS.out.tabs
        .combine(VARIANT_PREFILTER_TABLE.out.prefilter_table, by: 0)
        .combine(ch_priority_variants, by: 0)

    VAF_SPLIT_VARIANTS_HEXBIN(
        ch_vaf_split_input,
        vsRhoThr,
        vsLog10qThr,
        vsBinNv,
        vsBinNvShared,
        vsBinNr,
        vsBinVaf,
        vsBinVafAnchor,
        vsBinVafSingle,
        params.first_pass_betabinomial_cutoff,
        params.first_pass_binomial_cutoff
    )

    // -------------------------------------------------------------------------
    // Phase 14: Filter provenance — master variant × filter-status table
    // -------------------------------------------------------------------------
    ch_tab_nvnr_per_group = CREATE_TAB_NVNR.out
        .map { group, chr, mat_nv, mat_nr -> tuple(group, [mat_nv, mat_nr]) }
        .groupTuple(by: 0)
        .map { group, list -> tuple(group, list.flatten()) }

    ch_input_variant_provenance = FILTER_CHOSEN_VARIANTS_BY_BULK.out.bulk_filter_provenance
        .combine(CONCAT_FILTER_BINOM_BETABINOM_TAB_NV_NR.out.res_df,        by: 0)
        .combine(ch_vep_filter_provenance,                                  by: 0)
        .combine(ch_tab_nvnr_per_group,                                     by: 0)
        .combine(ch_usable_df_nv,                                           by: 0)
        .combine(ch_priority_variants,                                      by: 0)
        .combine(VARIANT_PREFILTER_TABLE.out.prefilter_table,               by: 0)
        .combine(VAF_SPLIT_VARIANTS_HEXBIN.out.all_outputs,                 by: 0)

    CUSTOM_VARIANT_FILTER_PROVENANCE(
        ch_input_variant_provenance
    )

    // -------------------------------------------------------------------------
    // Phase 15: Variant filter funnel (needs master_table + VAF cascade)
    // -------------------------------------------------------------------------
    ch_filter_funnel_input = CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table
        .combine(GET_VARIANTS_FROM_MERGED_VCF.out.all_variants, by: 0)
        .combine(VAF_SPLIT_VARIANTS_HEXBIN.out.cascade_binary,  by: 0)

    VARIANT_FILTER_FUNNEL(
        ch_filter_funnel_input
    )

    // -------------------------------------------------------------------------
    // Phase 16: Per-sample VCF splitting, focal pileup, and annotation
    // -------------------------------------------------------------------------
    LIST_SAMPLES_FROM_GROUP_VCF(
        ch_merged_vep_vcf_for_samples
    )

    // Focal pileup grid — fanned out per chromosome, then concatenated.
    ch_build_focal_input = CUSTOM_SOMATIC_SNPINDEL_FILTERRAWTABLES.out.pileup
        .combine(CUSTOM_VARIANT_FILTER_PROVENANCE.out.vcf_annotation_table, by: 0)
        .combine(LIST_SAMPLES_FROM_GROUP_VCF.out.sample_list,               by: 0)

    BUILD_FOCAL_PILEUP(
        ch_build_focal_input
    )

    CONCAT_FOCAL_PILEUP(
        BUILD_FOCAL_PILEUP.out.pileup_focal_chr.groupTuple(by: 0)
    )

    // One channel element per (group, sample_name).
    ch_per_sample = LIST_SAMPLES_FROM_GROUP_VCF.out.sample_list
        .map { group, f -> f.readLines().collect { s -> tuple(group, s) } }
        .flatMap { it }

    ch_annotate_input = ch_per_sample
        .combine(ch_merged_vep_vcf_for_samples,                             by: 0)
        .combine(CUSTOM_VARIANT_FILTER_PROVENANCE.out.vcf_annotation_table, by: 0)
        .combine(CONCAT_FOCAL_PILEUP.out.pileup_focal,                      by: 0)

    ANNOTATE_SAMPLE_VCF(
        ch_annotate_input
    )

    // -------------------------------------------------------------------------
    // Phase 17: Per-sample NR/NV/GT extraction → cross-sample matrices + schemes
    // -------------------------------------------------------------------------
    EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF(
        ANNOTATE_SAMPLE_VCF.out.annotated_vcf
    )

    ch_extracted_grouped = EXTRACT_NR_NV_GT_FROM_ANNOTATED_VCF.out.per_sample_vectors
        .groupTuple(by: 0)
        .map { group, sample_names, nr_files, nv_files, gt_files, vid_files ->
            tuple(group, sample_names, nr_files, nv_files, gt_files, vid_files)
        }

    ch_nr_nv_matrices_input = ch_extracted_grouped
        .combine(CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table, by: 0)

    CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF(
        ch_nr_nv_matrices_input,
        vsRhoThr,
        vsLog10qThr,
        vsBinNv,
        vsBinNvShared,
        vsBinNr,
        vsBinVaf,
        vsBinVafAnchor,
        vsBinVafSingle
    )

    // Tag each per-sample annotated VCF with scheme-membership flags, then merge
    // the per-sample scheme-annotated VCFs back into one group VCF.
    ch_scheme_annotate_input = ANNOTATE_SAMPLE_VCF.out.annotated_vcf
        .combine(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.scheme_membership, by: 0)

    ANNOTATE_VCF_SCHEME_MEMBERSHIP(
        ch_scheme_annotate_input
    )

    ch_merge_annotated_input = ANNOTATE_VCF_SCHEME_MEMBERSHIP.out.annotated_vcf
        .map { group, sample_name, vcf, tbi -> tuple(group, vcf, tbi) }
        .groupTuple(by: 0)

    MERGE_ANNOTATED_SAMPLE_VCFS(
        ch_merge_annotated_input
    )

    // Subset the merged scheme-annotated group VCF to the HQStat_QC_Depth scheme
    // (INFO/Scheme_HQStat_QC_Depth=1) — a compact, phylogeny-ready group VCF.
    SUBSET_ANNOTATED_VCF_HQSTAT_QC_DEPTH(
        MERGE_ANNOTATED_SAMPLE_VCFS.out.merged_vcf
    )

    // -------------------------------------------------------------------------
    // Phase 18: QC/filter diagnostics
    // -------------------------------------------------------------------------
    ch_hqstat_qc_input = CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table
        .join(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.scheme_membership, by: 0)

    COMPARE_HQSTAT_QC_REDUNDANCY(
        ch_hqstat_qc_input
    )

    ch_mandatory_qc_input = CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table
        .join(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.scheme_membership, by: 0)
        .join(MERGE_MANDATORY_PRIORITY_VARIANTS.out.merged_priority_variants, by: 0)

    MANDATORY_VARIANTS_QC_STATUS(
        ch_mandatory_qc_input
    )

    ch_downstream_collect = VAF_SPLIT_VARIANTS_HEXBIN.out.all_outputs
        .join(CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table,                by: 0)
        .join(MANDATORY_VARIANTS_QC_STATUS.out.status_table,                    by: 0)
        .join(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.scheme_membership,   by: 0)

    COLLECT_DOWNSTREAM_ARTIFACTS(
        ch_downstream_collect
    )

    // Lineage-inputs manifest for the downstream connector (one small task):
    //   index/lineage_inputs.csv       — "group,param,path" (multi-group safe)
    //   index/lineage_vcfs_<group>.csv — per-group VCF list (input_csv per group)
    // Keyed off phylogeny_matrices so only groups whose matrices were produced appear.
    def lineage_out_base = params.outputDir.toString().replaceAll('/+$', '') +
        "/workflow_outputs/${params.workspace}/${params.workflow_id}"

    ch_matrix_groups = VAF_SPLIT_VARIANTS_HEXBIN.out.phylogeny_matrices
        .map { group, _nv, _nr, _bin -> group }

    ch_manifest_in = ANNOTATE_SAMPLE_VCF.out.annotated_vcf
        .map { group, sample, _vcf, _tbi -> tuple(group, sample) }
        .groupTuple(by: 0)
        .join(ch_matrix_groups.map { g -> tuple(g, 'x') }, by: 0)
        .map { group, samples, _flag -> [group, samples] }
        .toList()

    EMIT_LINEAGE_INPUTS_MANIFEST( ch_manifest_in, lineage_out_base )

    ch_qc_threshold_input = ch_pileup_per_group
        .join(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.scheme_membership, by: 0)

    EXPLORE_QC_FILTER_THRESHOLDS(
        ch_qc_threshold_input,
        params.cutoff_as,
        params.cutoff_prop_clipped_reads,
        params.cutoff_prop_bp_under,
        params.cutoff_prop_bp_upper,
        params.cutoff_sd_indiv,
        params.cutoff_mad_indiv,
        params.cutoff_sd_both,
        params.cutoff_mad_both,
        params.cutoff_sd_extreme,
        params.cutoff_mad_extreme
    )

    ch_qc_influence_input = CONCAT_FOCAL_PILEUP.out.pileup_focal
        .join(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.scheme_membership, by: 0)

    QUANTIFY_QC_FILTER_INFLUENCE(
        ch_qc_influence_input
    )

    ch_cohort_metrics_input = VARIANT_FILTER_FUNNEL.out.filter_tracking
        .join(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.matrix_scheme_summary, by: 0)
        .join(COMPARE_HQSTAT_QC_REDUNDANCY.out.contingency,                       by: 0)
        .join(QUANTIFY_QC_FILTER_INFLUENCE.out.influence_tsv,                     by: 0)

    COHORT_METRICS_SUMMARY(
        ch_cohort_metrics_input
    )

    // -------------------------------------------------------------------------
    // Phase 19: Pileup-depth analysis + matrix scheme-summary plot
    // -------------------------------------------------------------------------
    ch_pileup_depth_mats = CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.nr_nv_matrices
        .map { group, nr_list, nv_list ->
            def schemes = params.nr_nv_pileup_depth_schemes.split(',').collect { it.trim() }.findAll { it }
            if (schemes.isEmpty()) {
                error "nr_nv_pileup_depth_schemes must list at least one scheme (group=${group})"
            }
            def nrs = nr_list instanceof List ? nr_list : [nr_list]
            def nvs = nv_list instanceof List ? nv_list : [nv_list]
            def nr_paths = []
            def nv_paths = []
            schemes.each { s ->
                def sfx = "_${s}.tsv"
                def nr_f = nrs.find { it.name.endsWith(sfx) }
                def nv_f = nvs.find { it.name.endsWith(sfx) }
                if (!nr_f || !nv_f) {
                    error "ANALYZE_NR_NV_PILEUP_DEPTH: missing NR or NV matrix for scheme '${s}' (group=${group})"
                }
                nr_paths.add(nr_f)
                nv_paths.add(nv_f)
            }
            tuple(group, nr_paths, nv_paths, schemes.join(','))
        }

    ANALYZE_NR_NV_PILEUP_DEPTH(
        ch_pileup_depth_mats,
        params.pileup_depth_min_sample_pct,
        params.pileup_depth_thresholds
    )

    PLOT_MATRIX_SCHEME_SUMMARY(
        CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.matrix_scheme_summary
            .join(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.matrix_per_sample_summary, by: 0)
            .join(CUSTOM_VARIANT_FILTER_PROVENANCE.out.upstream_per_sample, by: 0)
    )

    // =========================================================================
    // Phase 20: Germline and Allele Drop-Out (ADO) Analysis
    // =========================================================================
    ch_identify_germline_input = ch_usable_vcf_unpacked
        .combine(CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table, by: 0)
        .combine(ch_vep_filter_provenance,                          by: 0)
        .combine(FILTER_CHOSEN_VARIANTS_BY_BULK.out.bulk_filter_provenance, by: 0)

    IDENTIFY_GERMLINE_FROM_STATS(
        ch_identify_germline_input,
        params.germline_prev_pct
    )

    // Extract prevalence statistics per variant from annotated germline VCF
    EXTRACT_GERMLINE_PREVALENCE_TABLE(
        IDENTIFY_GERMLINE_FROM_STATS.out.annotated_vcf
    )

    // Generate distribution plots for germline prevalence
    PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS(
        EXTRACT_GERMLINE_PREVALENCE_TABLE.out.prevalence_table,
        params.germline_prev_pct
    )

    // Expand annotated group VCF to per-sample germline VCFs (3 filter sets each)
    // Re-use ch_per_sample (sample names already resolved by LIST_SAMPLES_FROM_GROUP_VCF)
    ch_subset_germline_input = ch_per_sample
        .combine(IDENTIFY_GERMLINE_FROM_STATS.out.annotated_vcf, by: 0)

    SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS(
        ch_subset_germline_input
    )

    // Create ADO tables from per-sample germline VCFs (3 parallel flavours)
    ch_ado_stats = SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS.out.stats_vcf
        .map { group, sample, vcf, tbi -> 
            def sample_ado_name = "${sample}_stats"
            tuple(sample_ado_name, vcf, tbi) 
        }

    ch_ado_vep = SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS.out.vep_vcf
        .map { group, sample, vcf, tbi -> 
            def sample_ado_name = "${sample}_vep"
            tuple(sample_ado_name, vcf, tbi) 
        }

    ch_ado_bulk = SUBSET_MERGED_VCF_HIGH_CONFIDENCE_GERMLINE_FROM_STATS.out.bulk_vcf
        .map { group, sample, vcf, tbi -> 
            def sample_ado_name = "${sample}_bulk"
            tuple(sample_ado_name, vcf, tbi) 
        }

    ch_all_ado_vcfs = ch_ado_stats.mix(ch_ado_vep, ch_ado_bulk)

    CREATE_ADO_TABLE_FROM_GERMLINE_VCF(
        ch_all_ado_vcfs,
        params.ado_sample_prop
    )

    // Summarize ADO intervals per sample — skip empty ADO tables (no variants in that filter set)
    SUMMARIZE_ADO_INTERVALS(
        CREATE_ADO_TABLE_FROM_GERMLINE_VCF.out.ado_table
            .filter { sample_name, tsv -> tsv.size() > 0 },
        params.ado_cov_cutoff
    )

    // Split res_ADO_* files by provenance label, then produce one labeled summary per filter set
    ch_df_sum_by_prov = SUMMARIZE_ADO_INTERVALS.out.df_sum
        .flatten()
        .branch {
            stats: it.name.contains('_stats')
            vep:   it.name.contains('_vep')
            bulk:  it.name.contains('_bulk')
        }

    CONCAT_ADO_STATS( ch_df_sum_by_prov.stats.collect(), 'stats' )
    CONCAT_ADO_VEP  ( ch_df_sum_by_prov.vep.collect(),   'vep'  )
    CONCAT_ADO_BULK ( ch_df_sum_by_prov.bulk.collect(),  'bulk' )

    ch_plot_ado_tables = CONCAT_ADO_STATS.out.merged_ADO
        .mix( CONCAT_ADO_VEP.out.merged_ADO  )
        .mix( CONCAT_ADO_BULK.out.merged_ADO )
        .collect()

    ch_plot_ado_summaries = CONCAT_ADO_STATS.out.summary_ADO
        .mix( CONCAT_ADO_VEP.out.summary_ADO  )
        .mix( CONCAT_ADO_BULK.out.summary_ADO )
        .collect()

    PLOT_ADO_GERMLINE_COMPARISON(
        ch_plot_ado_tables,
        ch_plot_ado_summaries
    )

    // Shared channel: (sample_name, [vcf, tbi]) for parquet + publish
    ch_annotated_vcf_per_sample = ANNOTATE_SAMPLE_VCF.out.annotated_vcf
        .map { group, sample_name, vcf, tbi -> tuple(sample_name, [vcf, tbi]) }

    // Convert somatic VCFs to Parquet for the lakehouse somatic_variants_summary table
    VCF_TO_PARQUET_SOMATIC_VARIANTS(
        ch_annotated_vcf_per_sample,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version
    )

    // Per-sample somatic variant counts from annotated VCFs → MultiQC generalstats TSVs
    ch_annotated_for_stats = ANNOTATE_SAMPLE_VCF.out.annotated_vcf
        .map { group, sample_name, vcf, tbi -> tuple(group, [vcf, tbi]) }
        .groupTuple(by: 0)
        .map { group, files -> tuple(group, files.flatten().collect()) }

    SOMATIC_STATS_FROM_ANNOTATED_VCF(
        ch_annotated_for_stats
    )

    // VEP per-chromosome summary HTMLs (ch_vep_html_for_multiqc defined in the VEP if/else above)

    ch_all_df_gt = PREPROCESS_VCF.out.df_gt
        .map { sample_name, df_gt, group -> df_gt }
        .collect()

    // -------------------------------------------------------------------------
    // Phase 12: Phylogenetic lineage analysis is handled by basej-lineage
    //           (Not part of somatic SNP/INDEL filtering pipeline)
    // -------------------------------------------------------------------------

    // Collect all MultiQC inputs and generate the somatic QC report
    MULTIQC_SOMATIC(
        SOMATIC_STATS_FROM_ANNOTATED_VCF.out.stats
            .flatten()
            .mix(SOMATIC_STATS_FROM_ANNOTATED_VCF.out.summary)
            .mix(ch_vep_html_for_multiqc)
            .mix(VARIANT_FILTER_FUNNEL.out.filter_plot_png
                .map { group, png -> png })
            .mix(PLOT_MATRIX_SCHEME_SUMMARY.out.scheme_summary_png
                .map { group, pngs -> pngs }.flatten())
            .mix(CREATE_NR_NV_MATRICES_FROM_ANNOTATED_VCF.out.first_round_hexbin
                .map { group, png -> png })
            .mix(VAF_SPLIT_VARIANTS_HEXBIN.out.all_outputs
                .map { group, files -> files }.flatten())
            .mix(COMPARE_HQSTAT_QC_REDUNDANCY.out.heatmap_png
                .map { group, png -> png })
            .mix(ANALYZE_NR_NV_PILEUP_DEPTH.out.sparsity_png
                .map { group, pngs -> pngs }.flatten())
            .mix(QUANTIFY_QC_FILTER_INFLUENCE.out.influence_png
                .map { group, png -> png })
            .mix(EXPLORE_QC_FILTER_THRESHOLDS.out.threshold_report
                .map { it -> it[2] }.flatten())
            .mix(COHORT_METRICS_SUMMARY.out.metrics_tsv
                .map { group, tsv -> tsv })
            .mix(EXTRACT_GERMLINE_PREVALENCE_TABLE.out.prevalence_table
                .map { group, tsv -> tsv })
            .mix(PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS.out.prevalence_plot
                .map { group, png -> png })
            .mix(PLOT_ADO_GERMLINE_COMPARISON.out.dist_plot)
            .mix(PLOT_ADO_GERMLINE_COMPARISON.out.summary_plot)
            .mix(PLOT_ADO_GERMLINE_COMPARISON.out.combined_plot)
            .mix(PLOT_ADO_GERMLINE_COMPARISON.out.summary_table)
            .collect(),
        params.project ?: 'N/A',
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        file("${projectDir}/assets/bioskryb_logo-tagline.png")
    )

    publish:
    // Primary outputs
    vcf_dnascope = SENTIEON_DNASCOPE.out.vcf
    
    vcf_somatic = ANNOTATE_SAMPLE_VCF.out.annotated_vcf
        .map { group, sample_name, vcf, tbi -> tuple(sample_name, vcf, tbi) }
    
    somatic_variants_summary = VCF_TO_PARQUET_SOMATIC_VARIANTS.out.parquet
    
    // Secondary analyses
    variant_tables = SEQUOIA_SECOND_FILTER_MERGE.out.df_filter
    
    provenance_report = CUSTOM_VARIANT_FILTER_PROVENANCE.out.master_table
    
    combined_report_pdf = VARIANT_FILTER_FUNNEL.out.filter_plot
        .map { group, pdf -> [groupId: group, pdf: pdf.toString()] }
    
    vep_annotated = ch_vep_vcf_for_publish
    
    germline_annotated_vcf = IDENTIFY_GERMLINE_FROM_STATS.out.annotated_vcf
    
    germline_prevalence_table = EXTRACT_GERMLINE_PREVALENCE_TABLE.out.prevalence_table
    
    germline_prevalence_plots = PLOT_GERMLINE_PREVALENCE_DISTRIBUTIONS.out.prevalence_plot
    
    ado_comparison_plots = PLOT_ADO_GERMLINE_COMPARISON.out.dist_plot
    
    // Filtering diagnostics + downstream handoff (June 2026 parity with Isai's branch)
    merged_annotated_vcf = MERGE_ANNOTATED_SAMPLE_VCFS.out.merged_vcf

    hqstat_qc_depth_vcf = SUBSET_ANNOTATED_VCF_HQSTAT_QC_DEPTH.out.subset_vcf

    downstream_artifacts = COLLECT_DOWNSTREAM_ARTIFACTS.out.bundle

    cohort_metrics = COHORT_METRICS_SUMMARY.out.metrics_tsv
        .mix(COHORT_METRICS_SUMMARY.out.metrics_json)

    filter_funnel = VARIANT_FILTER_FUNNEL.out.filter_tracking
        .mix(VARIANT_FILTER_FUNNEL.out.filter_report)
        .mix(VARIANT_FILTER_FUNNEL.out.filter_plot)

    // Phylogeny-ready matrices ("ForPhylogeny" scheme) — handoff to basej-lineage.
    // NR + NV feed SEQUOIA phylogeny/placement; the 0/1 binary matrix feeds the
    // mutational-signature stage (--binary_matrix). Emit order: NV, NR, binary.
    phylogeny_matrices = VAF_SPLIT_VARIANTS_HEXBIN.out.phylogeny_matrices
        .flatMap { group, nv, nr, bin -> [nv, nr, bin] }

    // Mandatory-variant QC status table — handoff to basej-lineage
    // (--mandatory_variants_qc_status); restricts heatmaps to passing mandatory variants.
    mandatory_qc_status = MANDATORY_VARIANTS_QC_STATUS.out.status_table
        .map { group, tsv -> tsv }

    lineage_inputs_manifest = EMIT_LINEAGE_INPUTS_MANIFEST.out.manifest

    lineage_vcf_csvs = EMIT_LINEAGE_INPUTS_MANIFEST.out.vcf_csvs.flatten()

    multiqc_report = MULTIQC_SOMATIC.out.report
}


// ============================================================================
// OUTPUT CONFIGURATION  (Seqera Platform)
// ============================================================================
output {

    // -------------------------------------------------------------------------
    // Primary analyses
    // -------------------------------------------------------------------------

    vcf_dnascope {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/variant_calls_dnascope"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "vcf_dnascope",
             tool:        "sentieon-dnascope",
             reference:   params.genome
    }

    vcf_somatic {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/vcf_somatic"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/vcf.csv"
            header(["biosampleName", "vcf", "tbi"])
        }
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "vcf_somatic",
             tool:        "basej-somatic",
             reference:   params.genome
    }

    // Lakehouse table — Hive-partitioned parquets under tables/somatic_variants_summary/
    somatic_variants_summary {
        path "tables"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "somatic_variants_summary",
             reference:   params.genome
    }

    // -------------------------------------------------------------------------
    // Secondary analyses
    // -------------------------------------------------------------------------

    variant_tables {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/variant_tables"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "variant_tables"
    }

    provenance_report {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/provenance"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "provenance_report",
             tool:        "custom-variant-filter-provenance"
    }

    combined_report_pdf {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/somatic_filter_reports"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/combined_report_pdf.csv"
            header true
        }
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "combined_report_pdf",
             tool:        "custom-variant-filter-provenance"
    }

    vep_annotated {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/vep_annotated"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "vep_annotated",
             tool:        "vep",
             reference:   params.genome
    }

    germline_annotated_vcf {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/germline_annotated"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "germline_annotated_vcf",
             tool:        "identify-germline-from-stats"
    }

    germline_prevalence_table {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/germline_prevalence"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "germline_prevalence_table"
    }

    germline_prevalence_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/germline_plots"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "germline_prevalence_plots"
    }

    ado_comparison_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/ado_plots"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "ado_comparison_plots",
             tool:        "ado-analysis"
    }

    // -------------------------------------------------------------------------
    // Filtering diagnostics + downstream handoff
    // -------------------------------------------------------------------------

    merged_annotated_vcf {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/merged_annotated_vcf"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "merged_annotated_vcf",
             reference:   params.genome
    }

    hqstat_qc_depth_vcf {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/hqstat_qc_depth_vcf"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "hqstat_qc_depth_vcf",
             reference:   params.genome
    }

    downstream_artifacts {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/downstream_inputs"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "downstream_artifacts"
    }

    cohort_metrics {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/cohort_metrics"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "cohort_metrics"
    }

    filter_funnel {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/filter_funnel"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "filter_funnel"
    }

    phylogeny_matrices {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/phylogeny_matrices"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "phylogeny_matrices"
    }

    mandatory_qc_status {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/mandatory_qc"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "mandatory_qc_status"
    }

    lineage_inputs_manifest {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/index"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "lineage_inputs_manifest"
    }

    lineage_vcf_csvs {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/index"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "lineage_vcf_csvs"
    }

    // -------------------------------------------------------------------------
    // Reports
    // -------------------------------------------------------------------------

    multiqc_report {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/reports"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "multiqc_report",
             tool:        "multiqc"
    }
}
