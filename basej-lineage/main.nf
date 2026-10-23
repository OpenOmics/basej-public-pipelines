nextflow.enable.dsl=2

// ============================================================================
// BASEJ-LINEAGE  —  matrix-consumer lineage + mutational-signatures pipeline
// ============================================================================
// Self-contained pipeline that MERGES two of Isai's standalone subworkflows:
//   * LINEAGE_FROM_MATRICES        — phylogeny / lineage from NR/NV matrices
//   * MUTSIGNATURES_FROM_MATRICES  — COSMIC mutational signatures from a binary matrix
//
// This thin orchestrator parses/validates params, builds the input channels,
// invokes both subworkflows, runs MULTIQC_LINEAGE, then publishes deliverables
// via the Seqera Platform `publish:` / `output {}` mechanism (no publishDir).
//
// Inputs (params):
//   nr_matrix / nv_matrix            (REQUIRED) — pre-built NR/NV PlusMandatoryNonEmpty TSVs
//                                                  (ForPhylogeny is placed in full by the
//                                                   phylogeny step → no heatmaps; use the
//                                                   fuller PlusMandatoryNonEmpty scheme)
//   binary_matrix                    (mutsig)   — 0/1 binary PlusMandatoryNonEmpty TSV
//   genotype_bin                     (optional) — pre-built binarized genotype matrix
//   mandatory_variants_qc_status     (optional) — mandatory-variant QC restriction
//   input_csv OR vcf_dir             (REQUIRED, one) — per-sample annotated VCFs
// ============================================================================

include { LINEAGE_FROM_MATRICES       } from './subworkflows/lineage_from_matrices.nf'
include { MUTSIGNATURES_FROM_MATRICES } from './subworkflows/mutsignatures_from_matrices.nf'

// ============================================================================
// MULTIQC — inline process (no publishDir; output handled via publish block)
// ============================================================================
process MULTIQC_LINEAGE {
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
custom_logo: ${logo.name}
custom_logo_width: 260

title: "basej-lineage v${pipeline_version}"
report_header_info:
  - Project: "${project}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false

run_modules:
  - custom_content

custom_content:
  order:
    - lineage_digital_heatmaps
    - lineage_signature_bargraphs
    - lineage_vaf_heatmaps

report_section_order:
  lineage_digital_heatmaps:
    order: 30
  lineage_signature_bargraphs:
    order: 20
  lineage_vaf_heatmaps:
    order: 10

EOF

    python3 << 'PYEOF'
import base64, glob

# Per Isai's requested lineage MultiQC layout — display exactly these figures, in order:
#   1) Digital genotype heatmaps: SNV, Indel, BOTH
#   2) Signature bar graphs:      SBS, DBS, ID
#   3) VAF heatmaps:              SNV, Indel, BOTH

def _img_b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

def _find(pattern):
    hits = glob.glob(pattern, recursive=True)
    return hits[0] if hits else None

def _write_png_section(fname, section_id, section_name, label_tpl, items):
    present = [(t, p) for t, p in items if p]
    if not present:
        return
    with open(fname, "w") as out:
        out.write(f"<!--\\nid: '{section_id}'\\nsection_name: '{section_name}'\\nplot_type: 'html'\\n-->\\n")
        out.write("<div style='padding:12px'>\\n")
        for t, p in present:
            out.write(f"<h3>{label_tpl.format(t=t)}</h3>\\n")
            out.write(f"<img src='data:image/png;base64,{_img_b64(p)}' style='width:100%; margin-bottom:16px; display:block;' />\\n")
        out.write("</div>\\n")

# ── 1) Digital genotype heatmaps (SNV, Indel, BOTH) ──
_write_png_section(
    "lineage_digital_heatmaps_mqc.html", "lineage_digital_heatmaps", "Digital Genotype Heatmaps",
    "Digital Genotype Heatmap ({t})",
    [(t, _find(f"**/res_composition_digital_{t}-1.png")) for t in ("SNV", "Indel", "BOTH")],
)

# ── 2) Signature bar graphs (SBS, DBS, ID) ──
_write_png_section(
    "lineage_signature_bargraphs_mqc.html", "lineage_signature_bargraphs", "Signature Bar Graphs",
    "Signature Bar Graph ({t})",
    [(t, _find(f"signature_bargraphs_{t}_combined.png")) for t in ("SBS", "DBS", "ID")],
)

# ── 3) VAF heatmaps (SNV, Indel, BOTH) ──
_write_png_section(
    "lineage_vaf_heatmaps_mqc.html", "lineage_vaf_heatmaps", "VAF Heatmaps",
    "VAF Heatmap ({t})",
    [(t, _find(f"**/res_composition_vaf_{t}-1.png")) for t in ("SNV", "Indel", "BOTH")],
)
PYEOF

    multiqc . \\
        --config multiqc_config.yaml \\
        --filename multiqc_report \\
        --force
    """
}

// ============================================================================
// MAIN WORKFLOW  (thin orchestrator)
// ============================================================================
workflow {
    main:

    // ── Parse / validate sequoia_phylogeny_mode ──────────────────────────────
    // Comma-separated subset of snv, indel, both, or "all".
    def raw_mode = params.sequoia_phylogeny_mode?.toString()?.trim()
    if ( !raw_mode ) {
        error "sequoia_phylogeny_mode cannot be empty"
    }
    def mode_parts = raw_mode.split(',').collect { it.trim().toLowerCase() }.findAll { it }
    if ( mode_parts.isEmpty() ) {
        error "sequoia_phylogeny_mode has no tokens after parsing: '${params.sequoia_phylogeny_mode}'"
    }
    def sequoia_phylo_modes
    if ( 'all' in mode_parts ) {
        if ( mode_parts.size() > 1 ) {
            log.warn "sequoia_phylogeny_mode: 'all' enables snv, indel, and both; extra tokens are ignored."
        }
        sequoia_phylo_modes = ['snv', 'indel', 'both'] as Set
    } else {
        def allowed = ['snv', 'indel', 'both']
        def bad = mode_parts.findAll { !(it in allowed) }
        if ( !bad.isEmpty() ) {
            error "sequoia_phylogeny_mode invalid token(s): ${bad.join(', ')}. Use snv, indel, both, and/or all (comma-separated)."
        }
        sequoia_phylo_modes = mode_parts as Set
    }
    def run_phylo_snv   = sequoia_phylo_modes.contains('snv')
    def run_phylo_indel = sequoia_phylo_modes.contains('indel')
    def run_phylo_both  = sequoia_phylo_modes.contains('both')

    // ── MPBoot bootstrap floor ────────────────────────────────────────────────
    // MPBoot hard-requires >= 1000 bootstrap replicates; a lower value silently
    // produces no .treefile and downstream read.tree() fails. Fail fast here.
    if ( (params.sequoia_mpboot_bootstrap as Integer) < 1000 ) {
        error "basej-lineage: sequoia_mpboot_bootstrap must be >= 1000 (MPBoot requirement); got ${params.sequoia_mpboot_bootstrap}."
    }

    // ── Lineage manifest parsing (--lineage_csv from somatic connector) ──────
    // When lineage_csv is provided, parse the CSV and populate params from its rows.
    // Format: group,param,path (single-group; first group wins, extra groups ignored).
    // Individual params (--nr_matrix, etc.) take precedence as explicit overrides.
    def lineage_csv_set = params.lineage_csv != null && !params.lineage_csv.toString().trim().isEmpty()
    def manifest_params = [:]  // param_name → path (from manifest CSV)
    if ( lineage_csv_set ) {
        log.info "basej-lineage: --lineage_csv provided — parsing manifest to resolve inputs."
        def csv_file = file(params.lineage_csv, checkIfExists: true)
        def lines = csv_file.text.trim().split('\n').toList()
        if ( lines.size() < 2 ) {
            error "basej-lineage: lineage_csv has no data rows (expected header: group,param,path)."
        }
        // Skip header, parse rows
        def target_group = null
        lines.drop(1).each { line ->
            def cols = line.split(',', 3)
            if ( cols.size() < 3 ) return  // skip malformed
            def g = cols[0].trim()
            def p = cols[1].trim()
            def v = cols[2].trim()
            // Single-group: lock to the first group encountered
            if ( target_group == null ) target_group = g
            if ( g != target_group ) {
                log.warn "basej-lineage: lineage_csv contains multiple groups; using '${target_group}', ignoring '${g}'."
                return
            }
            manifest_params[p] = v
        }
        log.info "basej-lineage: manifest resolved params for group '${target_group}': ${manifest_params.keySet().sort().join(', ')}"
    }

    // Helper: resolve a param — explicit CLI param wins, else manifest, else null/empty.
    def resolveParam = { String paramName, Object paramValue ->
        def explicitly_set = paramValue != null && !paramValue.toString().trim().isEmpty()
        if ( explicitly_set ) return paramValue.toString().trim()
        return manifest_params.containsKey(paramName) ? manifest_params[paramName] : null
    }

    // ── Required pre-built NR/NV matrices ─────────────────────────────────────
    def resolved_nr = resolveParam('nr_matrix', params.nr_matrix)
    def resolved_nv = resolveParam('nv_matrix', params.nv_matrix)
    if ( !resolved_nr || !resolved_nv ) {
        error "basej-lineage: both --nr_matrix and --nv_matrix are required (pre-built NR/NV matrices). Provide them directly or via --lineage_csv."
    }
    def nr_file = file(resolved_nr, checkIfExists: true)
    def nv_file = file(resolved_nv, checkIfExists: true)

    // ── Optional genotype_bin matrix (Option D) — empty file when absent ──
    def resolved_genotype_bin = resolveParam('genotype_bin', params.genotype_bin)
    def genotype_bin_set  = resolved_genotype_bin != null
    def genotype_bin_file = genotype_bin_set
        ? file(resolved_genotype_bin, checkIfExists: true)
        : file("${projectDir}/assets/dummy_file.txt")
    if ( genotype_bin_set ) {
        log.info "basej-lineage: genotype_bin provided (Option D) — phylogeny will use supplied genotype calls."
    } else {
        log.info "basej-lineage: no genotype_bin — phylogeny will discretize VAF from NR/NV."
    }

    // ── Optional mandatory_variants_qc_status — empty file when absent ──
    def resolved_mandatory_qc = resolveParam('mandatory_variants_qc_status', params.mandatory_variants_qc_status)
    def mandatory_qc_set  = resolved_mandatory_qc != null
    def mandatory_qc_file = mandatory_qc_set
        ? file(resolved_mandatory_qc, checkIfExists: true)
        : file("${projectDir}/assets/dummy_file.txt")
    if ( mandatory_qc_set ) {
        log.info "basej-lineage: mandatory_variants_qc_status provided — failing mandatory variants removed from heatmaps."
    }

    // ── VCF input (one of input_csv / vcf_dir) — for per-sample genotype tables ──
    def resolved_input_csv = resolveParam('input_csv', params.input_csv)
    def use_input_csv = resolved_input_csv != null
    def vcf_dir_set   = params.vcf_dir   != null && !params.vcf_dir.toString().trim().isEmpty()
    if ( use_input_csv && vcf_dir_set ) {
        log.warn "basej-lineage: both input_csv and vcf_dir are set; using input_csv."
    }
    if ( !use_input_csv && !vcf_dir_set ) {
        error "basej-lineage: set --input_csv (CSV with header: biosampleName,vcf [,vcf_index]) or --vcf_dir (dir of *_somatic_annotated.vcf.gz + .tbi). Provide them directly or via --lineage_csv. VCFs feed the per-sample genotype tables consumed by POSTPROCESS."
    }

    def ch_samples
    if ( use_input_csv ) {
        ch_samples = channel.fromPath( resolved_input_csv, checkIfExists: true )
            .splitCsv( header: true )
            .map { row ->
                def clean = row.collectEntries { k, v ->
                    [(k.toString().replace('\uFEFF', '').trim()): v?.toString()?.trim()]
                }
                def sample_name = clean.biosampleName
                def vcf_path    = clean.vcf
                if ( !sample_name || !vcf_path ) {
                    error "basej-lineage: each CSV row must define biosampleName and vcf (got keys: ${clean.keySet().sort().join(', ')})"
                }
                // Index is derived from the vcf path (.tbi) unless an explicit
                // vcf_index column is provided. Any extra columns (e.g. group) are ignored.
                def idx_path = ( clean.vcf_index && !clean.vcf_index.trim().isEmpty() )
                    ? clean.vcf_index
                    : "${vcf_path}.tbi"
                tuple( sample_name, file(vcf_path), file(idx_path) )
            }
            .ifEmpty { error "basej-lineage: no data rows in input_csv after the header" }
    } else {
        ch_samples = channel.fromFilePairs( "${params.vcf_dir}/*_somatic_annotated{,_filtered}.vcf.gz{,.tbi}", size: 2 )
            .map { _prefix, files ->
                def vcf = files.find { it.name.endsWith('.vcf.gz') && !it.name.endsWith('.tbi') }
                def tbi = files.find { it.name.endsWith('.tbi') }
                def sample_name = vcf.name.replaceAll('_somatic_annotated(_filtered)?\\.vcf\\.gz$', '')
                tuple( sample_name, vcf, tbi )
            }
            .ifEmpty { error "basej-lineage: no VCF+.tbi pairs matched in vcf_dir (expect *_somatic_annotated.vcf.gz)" }
    }

    // Per-sample channel keyed by cohort_id (single group): tuple(group, sample, vcf, tbi)
    def ch_per_sample = ch_samples.map { sample_name, vcf, tbi ->
        tuple( params.cohort_id, sample_name, vcf, tbi )
    }

    // ── Subworkflow 1: lineage / phylogeny from matrices ─────────────────────
    LINEAGE_FROM_MATRICES(
        ch_per_sample,
        params.cohort_id,
        nr_file,
        nv_file,
        genotype_bin_file,
        mandatory_qc_file,
        params.gender,
        run_phylo_snv,
        run_phylo_indel,
        run_phylo_both,
        params.sequoia_vaf_absent,
        params.sequoia_vaf_present,
        params.sequoia_tree_mut_pval,
        params.sequoia_keep_ancestral,
        params.sequoia_create_multi_tree,
        params.sequoia_genotype_conv_prob,
        params.sequoia_min_pval_for_true_somatic,
        params.sequoia_min_variant_reads_shared,
        params.sequoia_min_vaf_shared,
        params.sequoia_mpboot_path,
        params.sequoia_mpboot_bootstrap
    )

    // Fan the lineage_artifacts bundle: consumed by publish AND MultiQC.
    LINEAGE_FROM_MATRICES.out.lineage_artifacts
        .multiMap { it -> to_pub: it; to_mqc: it }
        .set { ch_lineage_bundle }

    // ── Subworkflow 2: mutational signatures from a binary matrix (optional) ──
    // Runs only when --binary_matrix is provided (or resolved from manifest);
    // otherwise the mutsig publish / MultiQC channels stay empty.
    def resolved_binary_matrix = resolveParam('binary_matrix', params.binary_matrix)
    def binary_set = resolved_binary_matrix != null

    ch_pub_sig_activities  = channel.empty()
    ch_pub_sig_bargraphs   = channel.empty()
    ch_pub_mutsig_bundle   = channel.empty()
    ch_pub_mutsig_coverage = channel.empty()
    ch_mqc_mutsig          = channel.empty()

    if ( binary_set ) {
        def ch_binary = channel.fromPath( resolved_binary_matrix, checkIfExists: true )

        MUTSIGNATURES_FROM_MATRICES(
            ch_binary,
            params.cohort_id,
            params.sig_genome_build,
            params.sig_context_type,
            params.sig_exome,
            params.sig_export_probabilities,
            params.sig_export_probabilities_per_mutation,
            params.musical_enabled,
            params.sigdyn_enabled,
            params.musicatk_enabled,
            params.musicatk_k_denovo
        )

        ch_pub_sig_activities = MUTSIGNATURES_FROM_MATRICES.out.activity_matrices.flatten()
            .mix(
                MUTSIGNATURES_FROM_MATRICES.out.musicatk_extra,
                MUTSIGNATURES_FROM_MATRICES.out.comparative_figs.flatten()
            )

        ch_pub_sig_bargraphs = MUTSIGNATURES_FROM_MATRICES.out.bargraph_figs.flatten()
            .mix( MUTSIGNATURES_FROM_MATRICES.out.bargraph_htmls.flatten() )

        ch_pub_mutsig_bundle   = MUTSIGNATURES_FROM_MATRICES.out.mutsig_artifacts
        ch_pub_mutsig_coverage = MUTSIGNATURES_FROM_MATRICES.out.coverage_tsvs.flatten()

        ch_mqc_mutsig = MUTSIGNATURES_FROM_MATRICES.out.bargraph_figs.flatten()
            .mix(
                MUTSIGNATURES_FROM_MATRICES.out.bargraph_htmls.flatten(),
                MUTSIGNATURES_FROM_MATRICES.out.coverage_tsvs.flatten(),
                MUTSIGNATURES_FROM_MATRICES.out.comparative_figs.flatten()
            )
    } else {
        log.info "basej-lineage: no --binary_matrix — mutational-signature branch skipped."
    }

    // ── MultiQC — aggregate key figures from both branches ───────────────────
    def ch_multiqc_inputs = channel.empty()
        .mix(
            ch_lineage_bundle.to_mqc.map { _group, bundle -> bundle },
            ch_mqc_mutsig
        )
        .collect()
        .ifEmpty( [] )

    MULTIQC_LINEAGE(
        ch_multiqc_inputs,
        params.project ?: 'N/A',
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        file("${projectDir}/assets/bioskryb_logo-tagline.png")
    )

    // =========================================================================
    // Publish channels — Seqera Platform output blocks
    // =========================================================================
    publish:

    // Genotype tables (per sample)
    genotype_tables = LINEAGE_FROM_MATRICES.out.genotype_table
        .map { group, sample_name, tsv -> [biosampleName: sample_name, group: group, genotype_table: tsv] }

    // Phylogeny outputs (per variant type)
    phylogeny_snv = LINEAGE_FROM_MATRICES.out.phylogeny_snv
        .flatMap { group, dirs -> (dirs instanceof List ? dirs : [dirs]).collect { [group: group, phylogeny: it] } }
    phylogeny_indel = LINEAGE_FROM_MATRICES.out.phylogeny_indel
        .flatMap { group, dirs -> (dirs instanceof List ? dirs : [dirs]).collect { [group: group, phylogeny: it] } }
    phylogeny_both = LINEAGE_FROM_MATRICES.out.phylogeny_both
        .flatMap { group, dirs -> (dirs instanceof List ? dirs : [dirs]).collect { [group: group, phylogeny: it] } }

    // Variant placement outputs (per variant type)
    placement_snv = LINEAGE_FROM_MATRICES.out.placement_snv
        .flatMap { group, dirs -> (dirs instanceof List ? dirs : [dirs]).collect { [group: group, placement: it] } }
    placement_indel = LINEAGE_FROM_MATRICES.out.placement_indel
        .flatMap { group, dirs -> (dirs instanceof List ? dirs : [dirs]).collect { [group: group, placement: it] } }
    placement_both = LINEAGE_FROM_MATRICES.out.placement_both
        .flatMap { group, dirs -> (dirs instanceof List ? dirs : [dirs]).collect { [group: group, placement: it] } }

    // Heatmaps (VAF + digital, per variant type) — POSTPROCESS emits individual files
    heatmaps = LINEAGE_FROM_MATRICES.out.heatmaps_snv
        .mix(LINEAGE_FROM_MATRICES.out.heatmaps_indel, LINEAGE_FROM_MATRICES.out.heatmaps_both)
        .flatMap { group, files -> (files instanceof List ? files : [files]).collect { [group: group, heatmaps: it] } }

    // Gathered lineage postprocess artifacts (figures + placed-variant tables + trees)
    lineage_artifacts = ch_lineage_bundle.to_pub
        .map { group, bundle -> [group: group, lineage_artifacts: bundle] }

    // ── Mutational signature (multi-tool) outputs — empty when mutsig skipped ──
    signature_activities = ch_pub_sig_activities
        .map { f -> [cohort: params.cohort_id, signature_activities: f] }
    signature_bargraphs = ch_pub_sig_bargraphs
        .map { f -> [cohort: params.cohort_id, signature_bargraph: f] }
    mutsig_artifacts = ch_pub_mutsig_bundle
        .map { cohort, bundle -> [cohort: cohort, mutsig_artifacts: bundle] }
    mutsig_coverage = ch_pub_mutsig_coverage
        .map { f -> [cohort: params.cohort_id, mutsig_coverage: f] }

    // MultiQC report
    multiqc_report = MULTIQC_LINEAGE.out.report
}


// ============================================================================
// OUTPUT CONFIGURATION  (Seqera Platform)
// ============================================================================
output {

    genotype_tables {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/genotype_tables"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/genotype_tables.csv"
            header true
        }
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "genotype_tables"
    }

    phylogeny_snv {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/phylogeny/snv"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "phylogeny_snv",
             tool:        "sequoia"
    }

    phylogeny_indel {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/phylogeny/indel"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "phylogeny_indel",
             tool:        "sequoia"
    }

    phylogeny_both {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/phylogeny/both"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "phylogeny_both",
             tool:        "sequoia"
    }

    placement_snv {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/variant_placement/snv"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "variant_placement_snv"
    }

    placement_indel {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/variant_placement/indel"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "variant_placement_indel"
    }

    placement_both {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/variant_placement/both"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "variant_placement_both"
    }

    heatmaps {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/heatmaps"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "heatmaps"
    }

    lineage_artifacts {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/lineage_artifacts"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "lineage_artifacts"
    }

    signature_activities {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/mutational_signatures/activities"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "signature_activities"
    }

    signature_bargraphs {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/mutational_signatures/bargraphs"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "signature_bargraphs"
    }

    mutsig_artifacts {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/mutational_signatures/artifacts"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "mutsig_artifacts"
    }

    mutsig_coverage {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/secondary_analyses/mutational_signatures/coverage"
        tags workspace:   params.workspace,
             project:     params.project,
             workflow_id: params.workflow_id,
             pipeline:    workflow.manifest.name,
             artifact:    "mutsig_coverage"
    }

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
