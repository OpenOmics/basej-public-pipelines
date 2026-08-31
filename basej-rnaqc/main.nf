nextflow.enable.dsl=2

// ============================================================================
// BASEJ-RNAQC: Single-Cell RNA QC Pipeline
// ============================================================================
// Description: FASTQ to STAR alignment, HTSeq quantification, QC metrics, and cluster QC
// Outputs: Per-biosample Parquet files for Athena queries + MultiQC report
// ============================================================================

// ============================================================================
// PROCESS: SAMTOOLS_SUBSAMPLE_CRAM
// Description: Calculate total reads and subsample CRAM to target (default 2M)
//              Uses same logic as SEQKIT_SAMPLE for consistency across platforms
// ============================================================================
process SAMTOOLS_SUBSAMPLE_CRAM {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(cram), path(crai), val(target_reads)
    val(samtools_seed)
    path(ref_fasta)

    output:
    tuple val(sample_name), path("${sample_name}.cram"), path("${sample_name}.cram.crai"), emit: cram
    tuple val(sample_name), env('TOTAL_READS'), env('FINAL_READS'), emit: read_counts
    path("${sample_name}_read_counts.txt"), emit: read_counts_file

    script:
    """
    # Count total reads in CRAM (reference required for CRAM decoding)
    export TOTAL_READS=\$(samtools view -c --reference '${ref_fasta}' '${cram}')
    echo "Total reads: \$TOTAL_READS"
    
    # Calculate subsample proportion (same logic as SEQKIT_SAMPLE)
    PROPORTION=\$(awk -v total="\$TOTAL_READS" -v target="${target_reads}" 'BEGIN { p=target/total; if(p>1) p=1; printf "%.8f", p }')
    echo "Subsample proportion: \$PROPORTION"
    
    # Construct seed argument: "seed.proportion"
    # Strip decimal from seed if present, and strip leading "0" from proportion
    SEED_RAW="${samtools_seed}"
    SEED_INT=\$(echo "\$SEED_RAW" | cut -d'.' -f1)
    PROP_FRAC=\$(echo "\$PROPORTION" | sed 's/^0//')
    SEED_ARG="\${SEED_INT}\${PROP_FRAC}"
    echo "SEED_RAW: \$SEED_RAW -> SEED_INT: \$SEED_INT, PROPORTION: \$PROPORTION -> PROP_FRAC: \$PROP_FRAC"
    echo "Seed argument: \$SEED_ARG"
    
    if [ \$(awk -v p="\$PROPORTION" 'BEGIN{if (p >= 1) print 1; else print 0}') -eq 1 ]; then
        echo "No subsampling needed..."
        # Check if input and output filenames are the same
        if [ '${cram}' != '${sample_name}.cram' ]; then
            ln -s '${cram}' '${sample_name}.cram'
        else
            echo "CRAM input already has correct name"
        fi
        if [ '${crai}' != '${sample_name}.cram.crai' ]; then
            ln -s '${crai}' '${sample_name}.cram.crai'
        else
            echo "CRAI input already has correct name"
        fi
        export FINAL_READS=\$TOTAL_READS
    else
        echo "Subsampling to ${target_reads} reads..."
        # Use samtools view -s for deterministic subsampling with seed
        # Write to temp file first to avoid reading/writing same file
        samtools view -s "\$SEED_ARG" --reference '${ref_fasta}' -C -o '${sample_name}_subsampled.cram' '${cram}'
        # If input has same name as output, remove originals before rename
        if [ '${cram}' = '${sample_name}.cram' ]; then
            rm -f '${sample_name}.cram' '${sample_name}.cram.crai'
        fi
        mv '${sample_name}_subsampled.cram' '${sample_name}.cram'
        samtools index '${sample_name}.cram' '${sample_name}.cram.crai'
        # Verify final read count
        export FINAL_READS=\$(samtools view -c --reference '${ref_fasta}' '${sample_name}.cram')
    fi
    
    echo "Final reads: \$FINAL_READS"
    
    # Write read counts to file for RNA_QC_PLOTS (matching SEQKIT_SAMPLE format)
    echo "\$TOTAL_READS" > '${sample_name}_read_counts.txt'
    echo "\$FINAL_READS" >> '${sample_name}_read_counts.txt'
    """
}

// ============================================================================
// PROCESS: MERGE_MULTILANE_FASTQ
// Description: Merge multi-lane FASTQs for the same biosample (cat R1s, cat R2s)
//              Only runs when read1/read2 contain pipe-delimited ("|") multi-lane paths
// ============================================================================
process MERGE_MULTILANE_FASTQ {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}_merged_R{1,2}.fastq.gz"), emit: reads

    script:
    """
    echo "Merging multi-lane FASTQs for sample: ${sample_name}"
    echo "R1 files: \$(ls *R1*.fastq.gz)"
    echo "R2 files: \$(ls *R2*.fastq.gz)"

    cat *R1*.fastq.gz > ${sample_name}_merged_R1.fastq.gz &
    cat *R2*.fastq.gz > ${sample_name}_merged_R2.fastq.gz
    wait

    echo "Merge complete for ${sample_name}"
    """
}

// ============================================================================
// PROCESS: SEQKIT_SAMPLE
// Description: Calculate total reads and subsample to target (default 2M)
//              Uses seqkit sample for ARM-native performance
// ============================================================================
process SEQKIT_SAMPLE {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(target_reads)
    val(seqkit_sample_seed)

    output:
    tuple val(sample_name), path("${sample_name}_subsampled_R*.fastq.gz"), emit: reads
    tuple val(sample_name), env('TOTAL_READS'), env('FINAL_READS'), emit: read_counts
    path("${sample_name}_read_counts.txt"), emit: read_counts_file

    script:
    def r1 = reads[0]
    def r2 = reads.size() == 2 ? reads[1] : ""
    def paired = reads.size() == 2
    def count_r2 = paired ? "( zcat '${r2}' | wc -l | awk '{print int(\$1/4)}' > read2.txt ) &" : ""
    def copy_r2 = paired ? "cp '${r2}' '${sample_name}_subsampled_R2.fastq.gz' &" : ""
    def sample_r2 = paired ? "seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j 2 -o '${sample_name}_subsampled_R2.fastq.gz' '${r2}'" : ""
    """
    set -euo pipefail

    # Count reads per mate in parallel
    ( zcat '${r1}' | wc -l | awk '{print int(\$1/4)}' > read1.txt ) &
    ${count_r2}
    wait

    R1=\$(cat read1.txt)
    if [ "${paired}" = "true" ]; then
      R2=\$(cat read2.txt)
      export TOTAL_READS=\$((R1 + R2))
    else
      export TOTAL_READS=\$R1
    fi
    echo "Total reads: \$TOTAL_READS"

    if [ "\$TOTAL_READS" -le "${target_reads}" ]; then
        echo "No subsampling needed, copying files..."
        cp '${r1}' '${sample_name}_subsampled_R1.fastq.gz' &
        ${copy_r2}
        wait
        export FINAL_READS=\$TOTAL_READS
    else
        echo "Subsampling to ${target_reads} reads with seqkit..."
        export PROPORTION=\$(awk -v t="${target_reads}" -v tot="\$TOTAL_READS" 'BEGIN { printf "%.18f", t/tot }')
        seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j 2 \\
          -o '${sample_name}_subsampled_R1.fastq.gz' '${r1}'
        ${paired ? sample_r2 : ""}
        export FINAL_READS=\$(zcat '${sample_name}_subsampled_R1.fastq.gz' | wc -l | awk '{print int(\$1/4)}')
    fi

    echo "Final reads: \$FINAL_READS"

    # Write read counts to file for RNA_QC_PLOTS
    echo "\$TOTAL_READS" > '${sample_name}_read_counts.txt'
    echo "\$FINAL_READS" >> '${sample_name}_read_counts.txt'
    """
}

// ============================================================================
// PROCESS: FASTP_TRIM
// Description: Adapter trimming with auto-detection
// ============================================================================
process FASTP_TRIM {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*_trim.fastq.gz"), emit: reads
    tuple val(sample_name), path("${sample_name}_fastp.json"), emit: json
    path("${sample_name}_fastp.json"), emit: json_flat

    script:
    """
    fastp --thread ${task.cpus} \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${sample_name}_R1_trim.fastq.gz --out2 ${sample_name}_R2_trim.fastq.gz \\
        --json ${sample_name}_fastp.json --html ${sample_name}_fastp.html \\
        --detect_adapter_for_pe \\
        2> ${sample_name}_fastp.log
    """
}

// ============================================================================
// PROCESS: STAR_ALIGN
// Description: STAR genome alignment with 2-pass mode
// ============================================================================
process STAR_ALIGN {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)
    path starindex

    output:
    tuple val(sample_name), path("star_outdir_${sample_name}/*Aligned.sortedByCoord.out.bam"), path("star_outdir_${sample_name}/*Chimeric.out.junction"), emit: bam
    path("star_outdir_${sample_name}/*Log.final.out"), emit: log_final
    path("star_outdir_${sample_name}"), emit: outdir

    script:
    def reads_paired_cmd = reads.size() == 2 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    """
    mkdir "star_outdir_${sample_name}"

    STAR \\
        --genomeDir ${starindex} \\
        --runThreadN ${task.cpus} \\
        --outSAMunmapped Within KeepPairs \\
        --twopassMode Basic \\
        --outReadsUnmapped None \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix "star_outdir_${sample_name}/${sample_name}_" \\
        --readFilesCommand zcat \\
        --readFilesIn ${reads_paired_cmd} \\
        --outSAMstrandField intronMotif \\
        --outSAMattributes All \\
        --outFilterMultimapNmax 10 \\
        --outSAMprimaryFlag OneBestScore \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 8 \\
        --chimOutJunctionFormat 1 \\
        --outSAMattrRGline "ID:${sample_name}"

    # Clean up intermediate files
    rm -rf star_outdir_${sample_name}/${sample_name}__STARgenome
    rm -rf star_outdir_${sample_name}/${sample_name}__STARpass1
    """
}

// ============================================================================
// PROCESS: SAMTOOLS_INDEX_FILTER
// Description: Index alignment and filter to primary alignments only
//              Handles both BAM from STAR and pre-aligned CRAM from Ultima
// ============================================================================
process SAMTOOLS_INDEX_FILTER {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(alignment), path(optional_file)
    path(ref_fasta)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam_bai

    script:
    def input_format = alignment.toString().endsWith('.cram') ? 'cram' : 'bam'
    def pl_tag       = input_format == 'cram' ? 'Ultima' : 'Illumina'
    """
    # Detect input format and process accordingly
    if [ "${input_format}" = "cram" ]; then
        echo "Input is CRAM format; converting to BAM for downstream compatibility"
        # Convert CRAM to BAM (reference required for CRAM decoding)
        # Primary-only filtering (-F 256) is applied in the shared pipe below
        samtools view -h -b --reference ${ref_fasta} -o ${sample_name}.raw.bam '${alignment}'
    else
        echo "Input is BAM format from STAR; filtering to primary alignments"

        # For STAR BAM, keep the input as-is (no index needed here; final BAM is indexed below)
        mv '${alignment}' ${sample_name}.raw.bam
    fi

    # Filter to primary alignments (-F 256) and update SM/PL tags in a single pass
    samtools view -F 256 -h ${sample_name}.raw.bam | \\
        sed -e 's/\\tSM:[^\\t]*//' -e 's/\\tPL:[^\\t]*//' | \\
        sed 's/^\\(@RG\\tID:[^\\t]*\\)/\\1\\tSM:${sample_name}\\tPL:${pl_tag}/' | \\
        samtools view -h -b - | \\
        samtools sort -O BAM -o ${sample_name}.bam
    
    # Index final BAM
    samtools index ${sample_name}.bam
    
    # Cleanup
    rm -f ${sample_name}.raw.bam
    """
}

// ============================================================================
// PROCESS: HTSEQ_COUNTS
// Description: HTSeq gene-level counting
// ============================================================================
process HTSEQ_COUNTS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(gtf_file)

    output:
    tuple val(sample_name), path("${sample_name}.htseq_counts.tsv"), emit: htseq_counts

    script:
    """
    htseq-count \\
        -f bam -r pos -s no -t exon -i gene_id --additional-attr=gene_name -m union \\
        --secondary-alignments=ignore \\
        ${bam} ${gtf_file} > ${sample_name}.htseq_counts.tsv
    """
}

// ============================================================================
// PROCESS: CREATE_HTSEQ_SUMMARY
// Description: Parse HTSeq counts to create gene summaries
// ============================================================================
process CREATE_HTSEQ_SUMMARY {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(htseq_file)
    path(tx2gene)

    output:
    path("df_gene_star_htseq_${sample_name}.tsv"), emit: gene_counts
    path("df_sum_detected_gene_${sample_name}.tsv"), emit: gene_types
    path("df_mtcounts_star_htseq_${sample_name}.tsv"), emit: mt_counts

    script:
    """
    echo -e "Working on ${htseq_file}"

    Rscript /usr/local/bin/htseq_to_gene_df.R ${htseq_file} ${tx2gene}

    # Rename outputs: R script creates files with .htseq_counts in the name
    # Need to remove .htseq_counts and ensure correct sample name
    for f in df_gene_star_htseq_*.tsv; do
        [ -e "\$f" ] || continue
        mv "\$f" df_gene_star_htseq_${sample_name}.tsv
    done

    for f in df_sum_detected_gene*.tsv; do
        [ -e "\$f" ] || continue
        mv "\$f" df_sum_detected_gene_${sample_name}.tsv
    done

    for f in df_mtcounts_star_htseq_*.tsv; do
        [ -e "\$f" ] || continue
        mv "\$f" df_mtcounts_star_htseq_${sample_name}.tsv
    done
    """
}

// ============================================================================
// PROCESS: MERGE_HTSEQ_SUMMARY
// Description: Merge all per-sample HTSeq summaries
// ============================================================================
process MERGE_HTSEQ_SUMMARY {
    tag "merge_htseq_summary"

    input:
    path(gene_counts_files)
    path(gene_types_files)
    path(mt_counts_files)

    output:
    path("df_gene_counts_starhtseq.tsv"), emit: gene_counts_tsv
    path("df_gene_types_detected_summary_starhtseq.tsv"), emit: gene_types_tsv
    path("df_mt_gene_counts_starhtseq.tsv"), emit: mt_counts_tsv
    tuple path("df_gene_counts_starhtseq.tsv"),
          path("df_gene_types_detected_summary_starhtseq.tsv"),
          path("df_mt_gene_counts_starhtseq.tsv"), emit: merge_tsv

    script:
    """
    # Merge gene counts
    echo -e "File\\tgene_id\\tcountHTSeq\\tgene_biotype\\tgene_symbol\\tgene_symbol_gene_id" > df_gene_counts_starhtseq.tsv
    cat df_gene_star_htseq_*.tsv | grep -v "^File" >> df_gene_counts_starhtseq.tsv

    # Merge gene types
    echo -e "File\\tgene_biotype\\tNumFeatures\\tPropFeatures\\tcountHTSeq\\tPropcountHTSeq" > df_gene_types_detected_summary_starhtseq.tsv
    cat df_sum_detected_gene_*.tsv | grep -v "^File" >> df_gene_types_detected_summary_starhtseq.tsv

    # Merge MT counts
    echo -e "File\\tTotalFeatures\\tMT_NumFeatures\\tMT_Counts\\tTotal_Counts\\tPropMT" > df_mt_gene_counts_starhtseq.tsv
    cat df_mtcounts_star_htseq_*.tsv | grep -v "^File" >> df_mt_gene_counts_starhtseq.tsv
    """
}

// ============================================================================
// PROCESS: CREATE_HTSEQ_MATRIX
// Description: Create expression matrix and housekeeping gene metrics,
//              then convert to long-format Parquet for lakehouse querying
// ============================================================================
process CREATE_HTSEQ_MATRIX {
    tag "create_htseq_matrix"

    input:
    tuple path(gene_counts), path(gene_types), path(mt_counts)
    val(workspace)
    val(workflow_id)
    val(dataset_id)
    val(pipeline_version)
    val(user)
    val(genome)

    output:
    path("matrix*"), emit: htseq_matrix
    path("HouseKeepingGenes_CV_mqc.tsv"), emit: housekeeping_genes_CV
    path("HouseKeepingGenes_Counts_mqc.tsv"), emit: housekeeping_genes_counts
    path("HKGenes_Expression__mqc.png"), emit: housekeeping_genes_clustergram
    path("gene_count_summary/workspace=*/workflow_id=*/output.parquet"), emit: gene_count_summary

    script:
    """
    Rscript /usr/local/bin/htseq_summarydf_to_matrix.R ${genome}

    python3 << 'PYEOF'
import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

workspace        = "${workspace}"
workflow_id      = "${workflow_id}"
dataset_id       = "${dataset_id}"
pipeline_version = "${pipeline_version}"
user             = "${user}"

# Read wide-format count matrix produced by htseq_summarydf_to_matrix.R
# Columns: gene_id_gene_symbol, <sample1>, <sample2>, ...
df_wide = pd.read_csv("matrix_gene_counts_starhtseq.txt", sep="\\t")

id_col      = "gene_id_gene_symbol"
sample_cols = [c for c in df_wide.columns if c != id_col]

# Melt to long (tidy) format: one row per (gene, sample)
df_long = df_wide.melt(
    id_vars=[id_col],
    value_vars=sample_cols,
    var_name="biosample",
    value_name="count"
)

# Split "ENSG00000001_GAPDH" → gene_id / gene_symbol on the last underscore
df_long["gene_id"]     = df_long[id_col].str.rsplit("_", n=1).str[0]
df_long["gene_symbol"] = df_long[id_col].str.rsplit("_", n=1).str[1]
df_long = df_long.drop(columns=[id_col])

# Attach run metadata
df_long["workspace"]        = workspace
df_long["workflow_id"]      = workflow_id
df_long["dataset_id"]       = dataset_id
df_long["pipeline"]         = "basej-rnaqc"
df_long["pipeline_version"] = pipeline_version
df_long["user"]             = user
df_long["count"]            = df_long["count"].astype(int)

df_long = df_long[[
    "gene_id", "gene_symbol", "biosample", "count",
    "workspace", "workflow_id", "dataset_id", "pipeline", "pipeline_version", "user"
]]

# Enforce schema types to match Iceberg gene_count_summary table
_str_cols    = ['gene_id','gene_symbol','biosample','workspace','workflow_id',
                'dataset_id','pipeline','pipeline_version','user']
_bigint_cols = ['count']
for _col in _bigint_cols:
    if _col in df_long.columns:
        df_long[_col] = pd.to_numeric(df_long[_col], errors='coerce').astype('Int64')
for _col in _str_cols:
    if _col in df_long.columns:
        df_long[_col] = df_long[_col].astype('string')

out_dir = f"gene_count_summary/workspace={workspace}/workflow_id={workflow_id}"
os.makedirs(out_dir, exist_ok=True)
pq.write_table(
    pa.Table.from_pandas(df_long, preserve_index=False),
    os.path.join(out_dir, "output.parquet")
)
print(f"Count matrix Parquet: {len(sample_cols)} samples x {len(df_wide)} genes → {len(df_long)} rows")
PYEOF
    """
}

// ============================================================================
// PROCESS: QUALIMAP_BAMRNA
// Description: Run Qualimap RNA-seq QC to analyze genomic origin of reads
// ============================================================================
process QUALIMAP_BAMRNA {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(gtf_file)

    output:
    path("qualimap_outdir_${sample_name}"), emit: outdir
    path("qualimap_outdir_${sample_name}/rnaseq_qc_results.txt"), emit: rnaseq_qc

    script:
    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp

    qualimap \\
        --java-mem-size=${task.memory.toGiga()}G \\
        rnaseq \\
        -gtf ${gtf_file} \\
        -bam ${bam} \\
        --sorted \\
        -outfile ${sample_name}.report.pdf \\
        -oc ${sample_name}.oc \\
        -outdir qualimap_outdir_${sample_name}
    """
}

// ============================================================================
// PROCESS: CREATE_QC_REPORT
// Description: Parse Qualimap output to extract genomic proportion metrics
// ============================================================================
process CREATE_QC_REPORT {
    tag "qc_report"

    input:
    path(qualimap_outdirs)

    output:
    path("qualimap_stats_mqc.csv"), emit: stats

    script:
    """
    echo "Working on Qualimap output directories"
    Rscript /usr/local/bin/parse_qualimap.R
    """
}

// ============================================================================
// PROCESS: PLOTTER_PCAHEATMAP
// Description: Generate cluster QC plots (PCA, Heatmap, MT%, GenesDetected)
// ============================================================================
process PLOTTER_PCAHEATMAP {
    tag "plotter_pcaheatmap"

    input:
    tuple path(gene_counts), path(gene_types), path(mt_counts)

    output:
    path("PCA_mqc.png"), emit: pca_plot
    path("Heatmap__mqc.png"), emit: heatmap_plot
    path("MT_mqc.png"), emit: mt_plot
    path("GenesDetected_mqc.png"), emit: genes_plot
    path("*.png"), emit: plots

    script:
    """
    Rscript /usr/local/bin/pca_heatmap_qc_salmon_htseq_nometadata.R

    # Rename files to uppercase for MultiQC compatibility
    mv pca_mqc.png PCA_mqc.png
    mv heatmap_mqc.png Heatmap__mqc.png
    mv mt_mqc.png MT_mqc.png
    mv genesdetected_mqc.png GenesDetected_mqc.png
    """
}

// ============================================================================
// PROCESS: RNA_QC_PLOTS
// Description: Consolidated mega-process - Parse metrics, create Parquet, generate plots
// ============================================================================
process RNA_QC_PLOTS {
    tag "rna_qc_plots"

    input:
    path(star_log_files)
    path(htseq_counts_files)
    path(fastp_jsons)
    path(read_counts_files)
    tuple path(gene_counts_tsv), path(gene_types_tsv), path(mt_counts_tsv)
    path(matrix_file)
    path(metadata_file)
    path(input_csv)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(user)

    output:
    path("rnaqc_summary/workspace=*/workflow_id=*/biosample=*/output.parquet"), emit: parquet
    path("*_selected_metrics_mqc.txt"), emit: mqc_metrics
    path("combined_selected_metrics.txt"), emit: combined_metrics
    path("rnaqc_all_metrics.tsv"), emit: summary_tsv
    path("composition_rnaqc.pdf"), emit: composition_pdf
    path("RNAQC_composition_mqc.jpg"), emit: composition_jpg
    path("RNA-QC_ConsensusScores.txt"), emit: consensus_scores
    path("RNA-QC_ConsensusScores_SummaryTable_mqc.txt"), emit: summary_verdict
    path("summary_verdict_group.txt"), emit: summary_verdict_group
    path("per_biosample_status.csv"), emit: per_biosample_status

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "RNA_QC_PLOTS: Starting consolidated QC workflow"
    echo "  1. Parse per-sample metrics → Parquets + TSVs"
    echo "  2. Merge all metrics into combined file"
    echo "  3. Generate visualizations"
    echo "=============================================="

    # ========== SECTION 1: Parse all metrics and create Parquets + TSVs ==========
    echo "Step 1: Parsing per-sample metrics and creating Parquets + TSVs..."

    python3 << 'PYEOF'
import os, json, glob, sys, re
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

dataset_id = "${dataset_id}"
workspace = "${workspace}"
workflow_id = "${workflow_id}"
pipeline_version = "${pipeline_version}"
user = "${user}"

# Discover samples from HTSeq counts files (works for both Illumina and Ultima)
# HTSeq counts are always produced regardless of platform
sample_names = set()
for f in glob.glob("*.htseq_counts.tsv"):
    # Extract sample name from: sample_name.htseq_counts.tsv
    sample = f.replace(".htseq_counts.tsv", "")
    if sample:
        sample_names.add(sample)

# Also discover from read_counts files (catches zero-read samples skipped in downstream)
for f in glob.glob("*_read_counts.txt"):
    sample = f.replace("_read_counts.txt", "")
    if sample:
        sample_names.add(sample)

# Fallback: try STAR log files (Illumina only)
if not sample_names:
    for f in glob.glob("*Log.final.out"):
        sample = f.replace("_Log.final.out", "").replace("Log.final.out", "")
        if sample:
            sample_names.add(sample)

if not sample_names:
    print("ERROR: No samples found (checked HTSeq counts and STAR logs)")
    sys.exit(1)

# Check if this is an Ultima run (no STAR logs)
has_star_logs = bool(glob.glob("*Log.final.out"))
platform = "Illumina" if has_star_logs else "Ultima"
print(f"Detected platform: {platform}")
print(f"Found {len(sample_names)} samples: {sorted(sample_names)}")

# Load merged HTSeq summary TSVs for gene detection metrics
df_gene_types = pd.read_csv("${gene_types_tsv}", sep="\\t")
df_mt_counts = pd.read_csv("${mt_counts_tsv}", sep="\\t")

all_summaries = []

# Process each sample
for sample in sorted(sample_names):
    summary = {
        "biosample": sample,
        "dataset_id": dataset_id,
        "pipeline": "basej-rnaqc",
        "pipeline_version": pipeline_version,
        "molecule_type": "rna",
        "user": user,
    }

    # ========== Parse Read Counts ==========
    read_count_files = glob.glob(f"{sample}*read_counts.txt")
    total_reads = 0
    final_reads = 0
    for f in read_count_files:
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                if len(lines) >= 2:
                    total_reads = int(lines[0].strip())
                    final_reads = int(lines[1].strip())
        except:
            pass

    summary["total_reads"] = total_reads
    summary["final_reads"] = final_reads
    summary["subsampled"] = total_reads != final_reads

    # ========== Parse Fastp JSON ==========
    fastp_files = glob.glob(f"{sample}*fastp.json")
    if fastp_files:
        try:
            with open(fastp_files[0], 'r') as f:
                fastp = json.load(f)
                before = fastp.get("summary", {}).get("before_filtering", {})
                summary["fastp_total_reads"] = before.get("total_reads", 0)
                summary["fastp_total_bases"] = before.get("total_bases", 0)
                summary["fastp_q20_rate"] = before.get("q20_rate", 0)
                summary["fastp_q30_rate"] = before.get("q30_rate", 0)
                summary["fastp_gc_content"] = before.get("gc_content", 0)
                summary["fastp_read1_mean_length"] = before.get("read1_mean_length", 0)
                summary["fastp_read2_mean_length"] = before.get("read2_mean_length", 0)
                after = fastp.get("summary", {}).get("after_filtering", {})
                summary["fastp_reads_after_filter"] = after.get("total_reads", 0)
                summary["fastp_q30_rate_after"] = after.get("q30_rate", 0)
                adapter = fastp.get("adapter_cutting", {})
                summary["fastp_adapter_trimmed_reads"] = adapter.get("adapter_trimmed_reads", 0)
                summary["fastp_adapter_trimmed_bases"] = adapter.get("adapter_trimmed_bases", 0)
                summary["fastp_duplication_rate"] = fastp.get("duplication", {}).get("rate", 0)
                filtering = fastp.get("filtering_result", {})
                total_before = before.get("total_reads", 0)
                passed = filtering.get("passed_filter_reads", 0)
                low_quality = filtering.get("low_quality_reads", 0)
                too_short = filtering.get("too_short_reads", 0)
                summary["fastp_pct_filtered"]    = (total_before - passed) / total_before if total_before > 0 else 0
                summary["fastp_pct_low_quality"] = low_quality / total_before if total_before > 0 else 0
                summary["fastp_pct_too_short"]   = too_short / total_before if total_before > 0 else 0
        except Exception as e:
            print(f"Warning: Could not parse fastp for {sample}: {e}")

    # ========== Parse STAR Log.final.out ==========
    star_log_files = glob.glob(f"*{sample}*Log.final.out")
    if star_log_files:
        try:
            with open(star_log_files[0], 'r') as f:
                star_data = {}
                for line in f:
                    if '|' in line:
                        parts = line.split('|')
                        if len(parts) == 2:
                            key = parts[0].strip()
                            value = parts[1].strip()
                            star_data[key] = value

                # Parse STAR metrics
                summary["star_input_reads"] = int(star_data.get("Number of input reads", "0"))
                summary["star_uniquely_mapped"] = int(star_data.get("Uniquely mapped reads number", "0"))
                summary["star_uniquely_mapped_pct"] = float(star_data.get("Uniquely mapped reads %", "0%").replace("%", ""))
                summary["star_avg_mapped_length"] = float(star_data.get("Average mapped length", "0"))
                summary["star_num_splices_total"] = int(star_data.get("Number of splices: Total", "0"))
                summary["star_num_splices_annotated"] = int(star_data.get("Number of splices: Annotated (sjdb)", "0"))
                summary["star_mismatch_rate"] = float(star_data.get("Mismatch rate per base, %", "0%").replace("%", ""))
                summary["star_deletion_rate"] = float(star_data.get("Deletion rate per base", "0%").replace("%", ""))
                summary["star_deletion_avg_length"] = float(star_data.get("Deletion average length", "0"))
                summary["star_insertion_rate"] = float(star_data.get("Insertion rate per base", "0%").replace("%", ""))
                summary["star_insertion_avg_length"] = float(star_data.get("Insertion average length", "0"))
                summary["star_multimapped"] = int(star_data.get("Number of reads mapped to multiple loci", "0"))
                summary["star_multimapped_pct"] = float(star_data.get("% of reads mapped to multiple loci", "0%").replace("%", ""))
                summary["star_multimapped_toomany"] = int(star_data.get("Number of reads mapped to too many loci", "0"))
                summary["star_unmapped_mismatches_pct"] = float(star_data.get("% of reads unmapped: too many mismatches", "0%").replace("%", ""))
                summary["star_unmapped_tooshort_pct"] = float(star_data.get("% of reads unmapped: too short", "0%").replace("%", ""))
                summary["star_unmapped_other_pct"] = float(star_data.get("% of reads unmapped: other", "0%").replace("%", ""))
                summary["star_chimeric_reads"] = int(star_data.get("Number of chimeric reads", "0"))
                summary["star_chimeric_pct"] = float(star_data.get("% of chimeric reads", "0%").replace("%", ""))
        except Exception as e:
            print(f"Warning: Could not parse STAR log for {sample}: {e}")

    # ========== Parse HTSeq Counts (assignment metrics from special rows) ==========
    htseq_files = glob.glob(f"{sample}*.htseq_counts.tsv")
    if htseq_files:
        try:
            htseq_assigned = 0
            htseq_no_feature = 0
            htseq_ambiguous = 0
            htseq_not_unique = 0
            htseq_total = 0

            with open(htseq_files[0], 'r') as f:
                for line in f:
                    parts = line.strip().split("\\t")
                    if len(parts) >= 2:
                        gene_id = parts[0]
                        count = int(parts[1]) if parts[1].isdigit() else 0

                        if gene_id.startswith("__"):
                            if gene_id == "__no_feature":
                                htseq_no_feature = count
                            elif gene_id == "__ambiguous":
                                htseq_ambiguous = count
                            elif gene_id == "__alignment_not_unique":
                                htseq_not_unique = count
                        else:
                            htseq_assigned += count

            htseq_total = htseq_assigned + htseq_no_feature + htseq_ambiguous + htseq_not_unique
            summary["htseq_total_reads"] = htseq_total
            summary["htseq_assigned"] = htseq_assigned
            summary["htseq_assigned_pct"] = (htseq_assigned / htseq_total) if htseq_total > 0 else 0
            summary["htseq_no_feature"] = htseq_no_feature
            summary["htseq_no_feature_pct"] = (htseq_no_feature / htseq_total) if htseq_total > 0 else 0
            summary["htseq_ambiguous"] = htseq_ambiguous
            summary["htseq_ambiguous_pct"] = (htseq_ambiguous / htseq_total) if htseq_total > 0 else 0
            summary["htseq_alignment_not_unique"] = htseq_not_unique
            summary["htseq_alignment_not_unique_pct"] = (htseq_not_unique / htseq_total) if htseq_total > 0 else 0
        except Exception as e:
            print(f"Warning: Could not parse HTSeq for {sample}: {e}")

    # ========== Parse Gene Detection Metrics from merged TSVs ==========
    try:
        # Gene types detected
        sample_gene_types = df_gene_types[df_gene_types["File"].str.contains(sample, na=False)]
        summary["genes_detected_total"] = int(sample_gene_types["NumFeatures"].sum()) if len(sample_gene_types) > 0 else 0
        summary["genes_protein_coding"] = int(sample_gene_types[sample_gene_types["gene_biotype"] == "protein_coding"]["NumFeatures"].iloc[0]) if len(sample_gene_types[sample_gene_types["gene_biotype"] == "protein_coding"]) > 0 else 0
        summary["genes_lncrna"] = int(sample_gene_types[sample_gene_types["gene_biotype"] == "lncRNA"]["NumFeatures"].iloc[0]) if len(sample_gene_types[sample_gene_types["gene_biotype"] == "lncRNA"]) > 0 else 0

        # MT metrics
        sample_mt = df_mt_counts[df_mt_counts["File"].str.contains(sample, na=False)]
        if len(sample_mt) > 0:
            summary["genes_mitochondrial"] = int(sample_mt["MT_NumFeatures"].iloc[0])
            summary["mt_gene_counts"] = int(sample_mt["MT_Counts"].iloc[0])
            summary["total_gene_counts"] = int(sample_mt["Total_Counts"].iloc[0])
            summary["mt_pct"] = float(sample_mt["PropMT"].iloc[0])
        else:
            summary["genes_mitochondrial"] = 0
            summary["mt_gene_counts"] = 0
            summary["total_gene_counts"] = 0
            summary["mt_pct"] = 0.0

        # rRNA (ribosomal) metrics — Qualimap does not emit a ribosomal metric, so we
        # derive it from the HTSeq gene-biotype summary (biotypes rRNA + Mt_rRNA).
        # Denominator = total assigned counts for the sample (same basis as mt_pct).
        rrna_biotypes = {"rRNA", "Mt_rRNA"}
        if len(sample_gene_types) > 0:
            rrna_rows = sample_gene_types[sample_gene_types["gene_biotype"].isin(rrna_biotypes)]
            rrna_counts = int(rrna_rows["countHTSeq"].sum())
            total_counts = int(sample_gene_types["countHTSeq"].sum())
            summary["genes_rrna"] = int(rrna_rows["NumFeatures"].sum())
            summary["rrna_gene_counts"] = rrna_counts
            summary["rrna_pct"] = (rrna_counts / total_counts) if total_counts > 0 else 0.0
        else:
            summary["genes_rrna"] = 0
            summary["rrna_gene_counts"] = 0
            summary["rrna_pct"] = 0.0
    except Exception as e:
        print(f"Warning: Could not parse gene detection for {sample}: {e}")

    all_summaries.append(summary)

    # ========== Save summary as JSON (parquets written after R scoring in Step 4) ==========
    with open(f"{sample}_summary.json", 'w') as jf:
        json.dump(summary, jf, default=str)

    # ========== Write individual TSV for MultiQC/R input (QC_Status filled in Step 4) ==========
    subsampled = "Yes" if summary.get('subsampled') else "No"

    with open(f"{sample}_selected_metrics_mqc.txt", 'w') as f:
        f.write("# id: 'rnaqc_summary'\\n")
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'QC Summary'\\n")
        f.write("# description: 'Per-sample alignment, expression, and genomic composition QC metrics with pass/fail status.'\\n")
        f.write("# pconfig:\\n")
        f.write("#   id: 'rnaqc_summary_table'\\n")
        f.write("sample_name\\tQC_Status\\tScore\\tMT_pct\\tProteinCodingGenes\\tTotalReads\\tFinalReads\\tfastp_q30_rate\\tPCT_Filtered\\tPCT_Low_Quality\\tPCT_Too_Short\\tPropIntronic\\n")
        f.write(f"{sample}\\tPENDING\\tNA\\t")
        f.write(f"{summary.get('mt_pct', 0):.4f}\\t{summary.get('genes_protein_coding', 0)}\\t")
        f.write(f"{total_reads}\\t{final_reads}\\t")
        f.write(f"{summary.get('fastp_q30_rate', 0):.4f}\\t")
        f.write(f"{summary.get('fastp_pct_filtered', 0):.4f}\\t")
        f.write(f"{summary.get('fastp_pct_low_quality', 0):.4f}\\t")
        f.write(f"{summary.get('fastp_pct_too_short', 0):.4f}\\t")
        f.write(f"NA\\n")

    print(f"  {sample}: metrics parsed, QC status pending R scoring")

print(f"Processed {len(all_summaries)} samples: created JSONs and TSVs")
PYEOF

    # ========== SECTION 2: Merge all per-sample TSV files ==========
    echo "Step 2: Merging per-sample metrics into combined file..."

    cat *_selected_metrics_mqc.txt | grep -v "^#" | grep -v "^sample_name" > temp_metrics.txt
    echo -e "sample_name\\tQC_Status\\tScore\\tMT%\\tProteinCodingGenes\\tTotalReads\\tFinalReads\\tfastp_q30_rate\\tPCT_Filtered\\tPCT_Low_Quality\\tPCT_Too_Short\\tPropIntronic" > combined_selected_metrics.txt
    cat temp_metrics.txt >> combined_selected_metrics.txt

    # Create CSV version for R script (rna_qc_plot.R expects comma-separated values)
    echo "sample_name,QC_Status,Score,MT%,ProteinCodingGenes,TotalReads,FinalReads,fastp_q30_rate,PCT_Filtered,PCT_Low_Quality,PCT_Too_Short,PropIntronic" > combined_selected_metrics.csv
    cat temp_metrics.txt | tr '\\t' ',' >> combined_selected_metrics.csv

    # ========== SECTION 3: Generate composition plots and consensus scores ==========
    echo "Step 3: Generating composition plots with Qualimap-derived genomic proportions..."

    # Parse Qualimap metadata to extract genomic proportions for each sample
    # The metadata_file is qualimap_stats_mqc.csv from CREATE_QC_REPORT
    # It contains RAW counts: "bam file", exonic, intronic, intergenic, reads aligned, not aligned, etc.

    python3 << 'PYEOF'
import pandas as pd
import re

# Read Qualimap metrics CSV (has raw counts, not percentages)
qualimap_df = pd.read_csv("${metadata_file}")

# Column name in CSV is "bam file" with a space, not "bam.file"
# Extract sample names from "bam file" column (remove path and .bam extension)
qualimap_df['biosampleName'] = qualimap_df['bam file'].apply(
    lambda x: re.sub(r'.*/', '', str(x)).replace('.bam', '').replace('_Aligned.sortedByCoord.out', '')
)

# Calculate genomic proportions from raw counts
# Total genomic features = exonic + intronic + intergenic + overlapping exon
qualimap_df['total_genomic'] = (
    qualimap_df['exonic'] +
    qualimap_df['intronic'] +
    qualimap_df['intergenic'] +
    qualimap_df['overlapping exon']
)

# Calculate proportions (0-1 range)
qualimap_df['prop_exonic'] = qualimap_df['exonic'] / qualimap_df['total_genomic']
qualimap_df['prop_intronic'] = qualimap_df['intronic'] / qualimap_df['total_genomic']
qualimap_df['prop_intergenic'] = qualimap_df['intergenic'] / qualimap_df['total_genomic']

# PropMappability = reads aligned / (reads aligned + not aligned)
qualimap_df['prop_mappability'] = (
    qualimap_df['reads aligned'] /
    (qualimap_df['reads aligned'] + qualimap_df['not aligned'])
)

# Read the combined metrics CSV to get MT% for each sample
metrics_df = pd.read_csv('combined_selected_metrics.csv')

# Merge MT% into qualimap_df based on sample name
qualimap_df = qualimap_df.merge(
    metrics_df[['sample_name', 'MT%']],
    left_on='biosampleName',
    right_on='sample_name',
    how='left'
)

# MT% is already a fraction (0-1) from the TSV
qualimap_df['prop_mt'] = qualimap_df['MT%']

# Read gene types TSV to get protein coding genes detected for each sample
gene_types_df = pd.read_csv('${gene_types_tsv}', sep='\\t')

# Filter for protein_coding biotype and extract NumFeatures (number of genes detected)
protein_coding_df = gene_types_df[gene_types_df['gene_biotype'] == 'protein_coding'].copy()

# Extract sample name from File column (remove .htseq_counts.tsv suffix)
protein_coding_df['biosampleName'] = protein_coding_df['File'].str.replace('.htseq_counts.tsv', '')

# Merge protein coding genes detected into qualimap_df
qualimap_df = qualimap_df.merge(
    protein_coding_df[['biosampleName', 'NumFeatures']],
    on='biosampleName',
    how='left'
)

# Rename NumFeatures to ProteinCodingGenesDetected for clarity
qualimap_df['ProteinCodingGenesDetected'] = qualimap_df['NumFeatures']

# Read input CSV to get groups information
input_csv_df = pd.read_csv('${input_csv}')

# Create a mapping of biosampleName to groups
sample_to_group = {}
if 'biosampleName' in input_csv_df.columns and 'groups' in input_csv_df.columns:
    sample_to_group = dict(zip(input_csv_df['biosampleName'], input_csv_df['groups']))

# Create enriched metadata with actual Qualimap values + MT% + Protein Coding Genes
metadata_rows = []
metadata_rows.append("SampleId,biosampleName,groups,PropExonic,PropIntergenic,PropMappability,PropIntronic,ProportionCountsMitochondrialGenes,ProteinCodingGenesDetected")

for _, row in qualimap_df.iterrows():
    sample = row['biosampleName']
    prop_exonic = row['prop_exonic']
    prop_intronic = row['prop_intronic']
    prop_intergenic = row['prop_intergenic']
    prop_mappability = row['prop_mappability']
    prop_mt = row['prop_mt']
    # A sample can have a Qualimap entry but no protein_coding row in the
    # gene-types summary, leaving NumFeatures as NaN after the left merge.
    # Fall back to 0 instead of crashing the whole process on int(NaN).
    _pcg = row['ProteinCodingGenesDetected']
    protein_coding_genes = int(_pcg) if pd.notna(_pcg) else 0

    # Get group from input CSV, default to 'C1' if not found
    group = sample_to_group.get(sample, 'C1')

    # Create metadata row with actual values
    metadata_rows.append(f"{sample},{sample},{group},{prop_exonic:.4f},{prop_intergenic:.4f},{prop_mappability:.4f},{prop_intronic:.4f},{prop_mt:.4f},{protein_coding_genes}")

# Write enriched metadata CSV
with open('metadata_enriched.csv', 'w') as f:
    f.write('\\n'.join(metadata_rows))

print(f"✓ Created metadata_enriched.csv with {len(metadata_rows)-1} samples using ACTUAL Qualimap-derived proportions + MT%")
print(f"  PropExonic range: {qualimap_df['prop_exonic'].min():.3f}-{qualimap_df['prop_exonic'].max():.3f}")
print(f"  PropIntronic range: {qualimap_df['prop_intronic'].min():.3f}-{qualimap_df['prop_intronic'].max():.3f}")
print(f"  PropIntergenic range: {qualimap_df['prop_intergenic'].min():.3f}-{qualimap_df['prop_intergenic'].max():.3f}")
print(f"  ProportionCountsMitochondrialGenes range: {qualimap_df['prop_mt'].min():.3f}-{qualimap_df['prop_mt'].max():.3f}")

# Write qualimap proportions for merging into parquet schema
qualimap_df[['biosampleName', 'prop_mappability', 'prop_exonic', 'prop_intergenic', 'prop_intronic']].rename(
    columns={'biosampleName': 'biosample'}
).to_csv('qualimap_proportions.tsv', sep='\t', index=False)
print(f"✓ Wrote qualimap_proportions.tsv with {len(qualimap_df)} samples")
PYEOF

    # Generate composition plots using rna_qc_plot.R with real Qualimap data
    echo "Calling rna_qc_plot.R with actual genomic proportion metrics..."
    Rscript /usr/local/bin/rna_qc_plot.R \\
        --matrix_file ${matrix_file} \\
        --metrics_file combined_selected_metrics.csv \\
        --metadata_file metadata_enriched.csv

    echo "✓ Composition plots generated successfully with Qualimap-derived metrics"

    echo ""
    echo "Step 4: Writing parquets and finalising QC status from R consensus scores..."
    python3 << 'PYEOF'
import os, json, glob
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

workspace   = "${workspace}"
workflow_id = "${workflow_id}"

# ========== Read R consensus scores ==========
df_scores = pd.read_csv("RNA-QC_ConsensusScores.txt", sep="\\t")
df_scores["CompositeScore"] = pd.to_numeric(df_scores["CompositeScore"], errors="coerce")
scores_by_sample = df_scores.set_index("SampleId")["CompositeScore"].to_dict()

def score_to_status(score):
    if score is None or pd.isna(score):
        return "FAIL"
    score = int(score)
    if score >= 4:
        return "PASS"
    elif score == 3:
        return "Borderline"
    else:
        return "FAIL"

# ========== Load qualimap proportions (scoring metrics) ==========
qualimap_props = {}
if os.path.exists('qualimap_proportions.tsv'):
    qdf = pd.read_csv('qualimap_proportions.tsv', sep='\t')
    for _, row in qdf.iterrows():
        qualimap_props[row['biosample']] = {
            'prop_mappability': row['prop_mappability'],
            'prop_exonic':      row['prop_exonic'],
            'prop_intergenic':  row['prop_intergenic'],
            'prop_intronic':    row['prop_intronic'],
        }

# ========== Schema column lists ==========
_str_cols    = ['biosample','dataset_id','pipeline','pipeline_version','molecule_type',
                'qc_status','workspace','workflow_id','user']
_double_cols = ['fastp_q20_rate','fastp_q30_rate','fastp_gc_content','fastp_q30_rate_after',
                'fastp_duplication_rate','star_uniquely_mapped_pct','star_avg_mapped_length',
                'star_mismatch_rate','star_deletion_rate','star_deletion_avg_length',
                'star_insertion_rate','star_insertion_avg_length','star_multimapped_pct',
                'star_unmapped_mismatches_pct','star_unmapped_tooshort_pct',
                'star_unmapped_other_pct','star_chimeric_pct','mt_pct',
                'htseq_assigned_pct','htseq_no_feature_pct','htseq_ambiguous_pct',
                'htseq_alignment_not_unique_pct',
                'prop_mappability','prop_exonic','prop_intergenic','prop_intronic','rrna_pct']
_bigint_cols = ['total_reads','final_reads','fastp_total_reads','fastp_total_bases',
                'fastp_read1_mean_length','fastp_read2_mean_length','fastp_reads_after_filter',
                'fastp_adapter_trimmed_reads','fastp_adapter_trimmed_bases',
                'star_input_reads','star_uniquely_mapped','star_num_splices_total',
                'star_num_splices_annotated','star_multimapped','star_multimapped_toomany',
                'star_chimeric_reads','htseq_total_reads','htseq_assigned',
                'htseq_no_feature','htseq_ambiguous','htseq_alignment_not_unique',
                'genes_detected_total','genes_protein_coding','genes_lncrna',
                'genes_mitochondrial','mt_gene_counts','total_gene_counts',
                'genes_rrna','rrna_gene_counts','qc_score']
_bool_cols   = ['subsampled']

# ========== Process each sample JSON ==========
json_files = sorted(glob.glob("*_summary.json"))
if not json_files:
    print("ERROR: No summary JSON files found")
    exit(1)

all_summaries = []
for jf in json_files:
    with open(jf) as fh:
        summary = json.load(fh)

    sample = summary["biosample"]
    composite_score = scores_by_sample.get(sample)
    qc_status = score_to_status(composite_score)

    # Override score for zero-read samples — ensure they show as FAIL / NA
    if int(summary.get("total_reads", 0) or 0) == 0:
        composite_score = None
        qc_status = "FAIL"

    summary["qc_status"] = qc_status
    summary["qc_score"]  = int(composite_score) if pd.notna(composite_score) else None

    # Merge qualimap proportions (scoring input metrics)
    if sample in qualimap_props:
        summary.update(qualimap_props[sample])
    else:
        summary.setdefault('prop_mappability', None)
        summary.setdefault('prop_exonic',      None)
        summary.setdefault('prop_intergenic',  None)
        summary.setdefault('prop_intronic',    None)

    # Write parquet
    df = pd.DataFrame([summary])
    for _col in _double_cols:
        if _col in df.columns:
            df[_col] = pd.to_numeric(df[_col], errors='coerce').astype('float64')
    for _col in _bigint_cols:
        if _col in df.columns:
            df[_col] = pd.to_numeric(df[_col], errors='coerce').astype('Int64')
    for _col in _bool_cols:
        if _col in df.columns:
            df[_col] = df[_col].astype('boolean')
    for _col in _str_cols:
        if _col in df.columns:
            df[_col] = df[_col].astype('string')
    out_dir = f"rnaqc_summary/workspace={workspace}/workflow_id={workflow_id}/biosample={sample}"
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   os.path.join(out_dir, "output.parquet"))

    all_summaries.append(summary)

    # Rewrite TSV with final QC_Status, Score, and all columns (including qualimap scoring metrics)
    score_str = f"{int(composite_score)}" if pd.notna(composite_score) else "NA"
    prop_map  = summary.get('prop_mappability')
    prop_exo  = summary.get('prop_exonic')
    prop_int  = summary.get('prop_intergenic')
    prop_intr = summary.get('prop_intronic')
    prop_map_str  = f"{prop_map:.4f}"  if prop_map  is not None else "NA"
    prop_exo_str  = f"{prop_exo:.4f}"  if prop_exo  is not None else "NA"
    prop_int_str  = f"{prop_int:.4f}"  if prop_int  is not None else "NA"
    prop_intr_str = f"{prop_intr:.4f}" if prop_intr is not None else "NA"
    with open(f"{sample}_selected_metrics_mqc.txt", 'w') as fh:
        fh.write("# id: 'rnaqc_summary'\\n")
        fh.write("# plot_type: 'table'\\n")
        fh.write("# section_name: 'QC Summary'\\n")
        fh.write("# description: 'Per-sample alignment, expression, and genomic composition QC metrics with pass/fail status.'\\n")
        fh.write("# pconfig:\\n")
        fh.write("#   id: 'rnaqc_summary_table'\\n")
        fh.write("sample_name\\tQC_Status\\tScore\\tPropMappability\\tPropExonic\\tPropIntergenic\\tPropIntronic\\tMT_pct\\tPropRibosomal\\tProteinCodingGenes\\tTotalReads\\tFinalReads\\tfastp_q30_rate\\tPCT_Filtered\\tPCT_Low_Quality\\tPCT_Too_Short\\n")
        fh.write(f"{sample}\\t{qc_status}\\t{score_str}\\t")
        fh.write(f"{prop_map_str}\\t{prop_exo_str}\\t{prop_int_str}\\t{prop_intr_str}\\t")
        fh.write(f"{summary.get('mt_pct', 0):.4f}\\t{summary.get('rrna_pct', 0):.4f}\\t{summary.get('genes_protein_coding', 0)}\\t")
        fh.write(f"{summary.get('total_reads', 0)}\\t{summary.get('final_reads', 0)}\\t")
        fh.write(f"{summary.get('fastp_q30_rate', 0):.4f}\\t")
        fh.write(f"{summary.get('fastp_pct_filtered', 0):.4f}\\t")
        fh.write(f"{summary.get('fastp_pct_low_quality', 0):.4f}\\t")
        fh.write(f"{summary.get('fastp_pct_too_short', 0):.4f}\\n")

    print(f"  {sample}: QC {qc_status} ({score_str}), CompositeScore={composite_score}")

# Write combined TSV of all samples (for nf-test validation)
pd.DataFrame(all_summaries).to_csv("rnaqc_all_metrics.tsv", sep="\\t", index=False)
pd.DataFrame(all_summaries)[['biosample','qc_status','pipeline','pipeline_version']].rename(
    columns={'biosample':'biosampleName'}).to_csv('per_biosample_status.csv', index=False)
print(f"Written {len(json_files)} parquets with R-derived QC status")
PYEOF
    echo "  ✓ Finished Step 4: parquets written with R consensus QC status"

    # Prepend MultiQC custom-content headers so summary_verdict.txt renders as a table in the report
    tmp=\$(mktemp)
    if [ -f summary_verdict.txt ] && [ ! -f RNA-QC_ConsensusScores_SummaryTable_mqc.txt ]; then
        mv summary_verdict.txt RNA-QC_ConsensusScores_SummaryTable_mqc.txt
    fi
    printf '# id: "qc_score_distribution"\n# section_name: "Total Usable Cells"\n# description: "Distribution of samples across composite QC score categories. Scores range from 0 (poor quality) to 5 (high quality)."\n# plot_type: "table"\n# pconfig:\n#   id: "qc_score_dist_table"\n#   title: "Total Usable Cells"\n' | cat - RNA-QC_ConsensusScores_SummaryTable_mqc.txt > "\$tmp" && mv "\$tmp" RNA-QC_ConsensusScores_SummaryTable_mqc.txt

    echo "RNA_QC_PLOTS: Successfully completed!"
    """
}

// ============================================================================
// PROCESS: MULTIQC
// Description: Aggregate all metrics and plots into final HTML report
// ============================================================================
process MULTIQC {
    tag "multiqc"

    input:
    path(mqc_files)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    path(logo)

    output:
    path("multiqc_report.html"), emit: report

    script:
    """
    cat > multiqc_config_runtime.yaml << EOF
custom_logo_title: 'BioSkryb Genomics'
custom_logo: bioskryb_logo-tagline.png
custom_logo_width: 260

title: "basej-rnaqc v${pipeline_version}"
report_header_info:
  - Dataset ID: "${dataset_id}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false
skip_generalstats: true

module_order:
  - custom_content
  - fastp

report_section_order:
  rnaqc_summary:
    order: 1000
  RNAQC_composition:
    order: 900
  qc_score_distribution:
    order: 800

fn_clean_exts:
  - ".csv"
  - ".tsv"
  - "_mqc"

table_cond_formatting_rules:
  QC_Status:
    pass:
      - s_eq: "PASS"
    warn:
      - s_eq: "Borderline"
    fail:
      - s_eq: "FAIL"

custom_data:
  Heatmap:
    section_name: "Heatmap "
  RNAQC_composition:
    section_name: "QC Composition"
    description: "Per-sample QC metric distributions grouped by cluster. Each panel shows a key metric across all cells."
  rnaqc_summary:
    headers:
      QC_Status:
        title: "QC Status"
        description: "Composite QC status: PASS, Borderline, or FAIL"
        placement: 100
      Score:
        title: "Score"
        description: "Composite QC score (0-5, e.g. 3)"
        placement: 110
      PropMappability:
        title: "Proportion Mappable"
        description: "Proportion of reads mapping to the genome (0-1). Source: Qualimap."
        format: "{:.4f}"
        placement: 120
      PropExonic:
        title: "Proportion Exonic"
        description: "Proportion of mapped reads falling in exonic regions (0-1). Source: Qualimap."
        format: "{:.4f}"
        placement: 130
      PropIntergenic:
        title: "Proportion Intergenic"
        description: "Proportion of mapped reads falling in intergenic regions (0-1). Source: Qualimap."
        format: "{:.4f}"
        placement: 140
      PropIntronic:
        title: "Proportion Intronic"
        description: "Proportion of mapped reads falling in intronic regions (0-1). Source: Qualimap."
        format: "{:.4f}"
        placement: 195
      MT_pct:
        title: "Proportion Mitochondrial"
        description: "Proportion of reads mapping to the mitochondrial genome (0-1)."
        format: "{:.4f}"
        placement: 150
      PropRibosomal:
        title: "Proportion Ribosomal"
        description: "Proportion of assigned counts from ribosomal RNA genes (rRNA + Mt_rRNA biotypes) (0-1). Source: HTSeq gene-biotype counts."
        format: "{:.4f}"
        placement: 155
      ProteinCodingGenes:
        title: "Protein Coding Genes"
        description: "Number of protein-coding genes detected (count > 0). Source: DESeq2/featureCounts."
        format: "{:,.0f}"
        placement: 160
      TotalReads:
        title: "Total Reads"
        description: "Total read pairs before trimming and filtering."
        format: "{:,.0f}"
        placement: 170
      FinalReads:
        title: "Final Reads"
        description: "Read pairs after fastp trimming."
        format: "{:,.0f}"
        placement: 180
      fastp_q30_rate:
        title: "Proportion Q30"
        description: "Proportion of bases with Phred quality >= 30 after fastp trimming (0-1). Source: fastp."
        format: "{:.4f}"
        placement: 190
      PCT_Filtered:
        title: "Proportion Filtered"
        description: "Proportion of reads removed by fastp (all filters combined: low quality, too short, too many N) (0-1). Source: fastp filtering_result."
        format: "{:.4f}"
        placement: 200
      PCT_Low_Quality:
        title: "Proportion Low Quality"
        description: "Proportion of reads removed due to low base quality (0-1). Source: fastp filtering_result.low_quality_reads."
        format: "{:.4f}"
        placement: 210
      PCT_Too_Short:
        title: "Proportion Too Short"
        description: "Proportion of reads removed for being too short after trimming (0-1). Source: fastp filtering_result.too_short_reads."
        format: "{:.4f}"
        placement: 220
EOF

    multiqc . -n multiqc_report.html -c multiqc_config_runtime.yaml --force
    """
}

// ============================================================================
// WORKFLOW ORCHESTRATION
// ============================================================================
workflow {
    main:
    // Validate required parameters
    if (!params.input_csv) {
        exit 1, "ERROR: --input_csv parameter is required"
    }

    if (!params.gtf) {
        exit 1, "ERROR: GTF file not found for genome ${params.genome}. Check genomes.config"
    }

    if (!params.tx2gene) {
        exit 1, "ERROR: tx2gene file not found for genome ${params.genome}. Check genomes.config"
    }

    // Auto-detect platform from CSV: parse with splitCsv and check whether
    // *any* row has a non-empty 'cram' column (handles blank first rows and
    // quoted fields that would trip up a naive split(',') approach).
    def csv_file = file(params.input_csv)
    def csv_rows = csv_file.splitCsv(header: true)
    def is_ultima = csv_rows.any { row -> row.cram?.trim() }

    log.info "Auto-detected platform: ${is_ultima ? 'Ultima (CRAM)' : 'Illumina (FASTQ)'}"

    if (is_ultima) {
        // ========== Ultima path: CRAM input (pre-aligned) ==========
        ch_reads = channel.fromPath(params.input_csv, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row.cram }
            .map { row ->
                def cram = file(row.cram)
                def crai = row.crai ? file(row.crai) : file(row.cram + '.crai')
                [row.biosampleName, cram, crai]
            }

        ch_reads.view { sample -> "Processing Ultima sample (pre-aligned CRAM): ${sample[0]}" }
        ch_reads.ifEmpty { exit 1, "ERROR: No CRAM files specified in --input_csv" }

        // CRAM subsampling
        SAMTOOLS_SUBSAMPLE_CRAM(
            ch_reads.map { sample_id, cram_file, crai_file -> [sample_id, cram_file, crai_file, params.n_reads] },
            params.samtools_seed,
            file(params.ref_fasta, checkIfExists: true)
        )

        // Pass subsampled CRAMs to INDEX_FILTER
        ch_for_index = SAMTOOLS_SUBSAMPLE_CRAM.out.cram
            .map { sample_id, cram_file, crai_file -> [sample_id, cram_file, crai_file] }

        // Store read count metrics for RNA_QC_PLOTS
        ch_readcount_metrics = SAMTOOLS_SUBSAMPLE_CRAM.out.read_counts_file.collect()
        // Create empty fastp metrics (not applicable for Ultima)
        ch_fastp_metrics = Channel.value([])
        // Create empty STAR log files (not applicable for Ultima)
        ch_star_logs = Channel.value([])

        // Per-sample read counts for zero-read filtering
        ch_read_counts = SAMTOOLS_SUBSAMPLE_CRAM.out.read_counts

    } else {
        // ========== Illumina path: FASTQ input ==========
        // Multi-lane support: read1/read2 can contain pipe-delimited ("|") paths
        //   Single-lane: biosampleName,s3://bucket/R1.fastq.gz,s3://bucket/R2.fastq.gz
        //   Multi-lane:  biosampleName,s3://bucket/L001_R1.fq.gz|s3://bucket/L002_R1.fq.gz,s3://bucket/L001_R2.fq.gz|s3://bucket/L002_R2.fq.gz
        ch_reads_csv = channel.fromPath(params.input_csv, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row.read1 && row.read2 }

        // Branch on pipe character to detect multi-lane vs single-lane
        ch_reads_branched = ch_reads_csv
            .branch {
                multilane: it.read1.contains('|')
                singlelane: true
            }

    // Single-lane: one R1, one R2 — pass through directly
    ch_singlelane = ch_reads_branched.singlelane
        .map { row -> [row.biosampleName, [file(row.read1, checkIfExists: true), file(row.read2, checkIfExists: true)]] }

    ch_singlelane.view { sample -> "Processing Illumina sample: ${sample[0]}" }

    // Multi-lane: split on "|", collect file objects, flatten for cat
    MERGE_MULTILANE_FASTQ(
        ch_reads_branched.multilane
            .map { row ->
                def r1_files = row.read1.tokenize('|').collect { file(it.trim(), checkIfExists: true) }
                def r2_files = row.read2.tokenize('|').collect { file(it.trim(), checkIfExists: true) }
                [row.biosampleName, r1_files + r2_files]
            }
    )

    ch_fastq_input = ch_singlelane.mix(MERGE_MULTILANE_FASTQ.out.reads)

        // FASTQ preprocessing pipeline
        SEQKIT_SAMPLE(
        ch_fastq_input.map { sample_id, reads -> [sample_id, reads, params.n_reads] 
        },
        params.seqkit_sample_seed
    )
        FASTP_TRIM(SEQKIT_SAMPLE.out.reads)

        // Alignment
        if (!params.star_index) {
            exit 1, "ERROR: STAR index not found for genome ${params.genome}. Check genomes.config"
        }
        STAR_ALIGN(FASTP_TRIM.out.reads, file(params.star_index, checkIfExists: true))

        // Convert STAR BAM output to channel format for INDEX_FILTER
        ch_for_index = STAR_ALIGN.out.bam
            .map { sample_id, bam_path, junction_path -> [sample_id, bam_path, junction_path] }

        // Store metrics
        ch_readcount_metrics = SEQKIT_SAMPLE.out.read_counts_file.collect()
        ch_fastp_metrics = FASTP_TRIM.out.json_flat.collect()
        ch_star_logs = STAR_ALIGN.out.log_final.collect()

        // Per-sample read counts for zero-read filtering
        ch_read_counts = SEQKIT_SAMPLE.out.read_counts
    }

    // ========== Common downstream processing (both platforms) ==========
    // Filter out zero-read samples before heavy processing
    // Zero-read samples will still appear in RNA_QC_PLOTS with 0/NA values (discovered via read_counts files)
    ch_for_index_with_counts = ch_for_index
        .join(ch_read_counts)   // [sample, alignment, optional_file, total_reads, final_reads]

    ch_for_index_nonzero = ch_for_index_with_counts
        .filter { sample, alignment, optional_file, total_reads, final_reads ->
            if (total_reads.toString().toLong() == 0) {
                log.warn "Skipping downstream processing for sample '${sample}' — 0 reads"
                return false
            }
            return true
        }
        .map { sample, alignment, optional_file, total_reads, final_reads -> [sample, alignment, optional_file] }

    // Index, filter, and convert CRAM to BAM (Ultima) or filter BAM (Illumina)
    SAMTOOLS_INDEX_FILTER(ch_for_index_nonzero, file(params.ref_fasta, checkIfExists: true))

    // Quantification and QC
    QUALIMAP_BAMRNA(SAMTOOLS_INDEX_FILTER.out.bam_bai, file(params.gtf, checkIfExists: true))
    HTSEQ_COUNTS(SAMTOOLS_INDEX_FILTER.out.bam_bai, file(params.gtf, checkIfExists: true))

    // ========== HTSeq Processing ==========
    CREATE_HTSEQ_SUMMARY(
        HTSEQ_COUNTS.out.htseq_counts,
        file(params.tx2gene, checkIfExists: true)
    )

    MERGE_HTSEQ_SUMMARY(
        CREATE_HTSEQ_SUMMARY.out.gene_counts.collect(),
        CREATE_HTSEQ_SUMMARY.out.gene_types.collect(),
        CREATE_HTSEQ_SUMMARY.out.mt_counts.collect()
    )

    CREATE_HTSEQ_MATRIX(
        MERGE_HTSEQ_SUMMARY.out.merge_tsv,
        params.workspace,
        params.workflow_id,
        params.dataset_id,
        workflow.manifest.version,
        params.pipeline_user,
        params.genome
    )

    // ========== QC Analysis ==========
    PLOTTER_PCAHEATMAP(MERGE_HTSEQ_SUMMARY.out.merge_tsv)

    // Parse Qualimap outputs to extract genomic proportion metrics
    CREATE_QC_REPORT(QUALIMAP_BAMRNA.out.outdir.collect())

    // ========== Consolidated RNA QC + Parquet Generation ==========
    RNA_QC_PLOTS(
        ch_star_logs,
        HTSEQ_COUNTS.out.htseq_counts.map { it[1] }.collect(),  // Extract just the file path from tuple
        ch_fastp_metrics,
        ch_readcount_metrics,
        MERGE_HTSEQ_SUMMARY.out.merge_tsv,
        CREATE_HTSEQ_MATRIX.out.htseq_matrix,
        CREATE_QC_REPORT.out.stats,
        file(params.input_csv, checkIfExists: true),
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        params.pipeline_user
    )

    // ========== MultiQC Report ==========
    // NOTE: mix() already tolerates empty channels, so optional/ignored-process
    // outputs are mixed directly. Do NOT use .ifEmpty(channel.empty()) here:
    // ifEmpty expects a *value*, and passing a channel makes an empty input
    // emit the channel object itself, which then reaches MULTIQC as a bogus
    // path ("Not a valid path value type: ... DataflowStream").
    ch_mqc = RNA_QC_PLOTS.out.mqc_metrics.flatten()
        .mix(ch_fastp_metrics.flatten())
        .mix(PLOTTER_PCAHEATMAP.out.plots.flatten())
        .mix(RNA_QC_PLOTS.out.composition_jpg)
        .mix(CREATE_HTSEQ_MATRIX.out.housekeeping_genes_CV)
        .mix(CREATE_HTSEQ_MATRIX.out.housekeeping_genes_counts)
        .mix(CREATE_HTSEQ_MATRIX.out.housekeeping_genes_clustergram)
        .mix(RNA_QC_PLOTS.out.summary_verdict)
        .collect()

    MULTIQC(
        ch_mqc,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        file("${projectDir}/assets/bioskryb_logo-tagline.png")
    )

    publish:
    bam_files = SAMTOOLS_INDEX_FILTER.out.bam_bai
        .map { sample_name, bam, bai -> [biosampleName: sample_name, bam: bam, bai: bai] }

    pca_heatmap_plots = channel.empty()
        .mix(PLOTTER_PCAHEATMAP.out.pca_plot)
        .mix(PLOTTER_PCAHEATMAP.out.heatmap_plot)
        .mix(PLOTTER_PCAHEATMAP.out.mt_plot)
        .mix(PLOTTER_PCAHEATMAP.out.genes_plot)
        .collect()
        .ifEmpty([])

    rnaqc_summary = RNA_QC_PLOTS.out.parquet

    gene_count_summary = CREATE_HTSEQ_MATRIX.out.gene_count_summary

    rnaqc_metrics_combined = RNA_QC_PLOTS.out.combined_metrics

    // Full-schema TSV of all samples for nf-test validation
    rnaqc_summary_tsv = RNA_QC_PLOTS.out.summary_tsv

    rnaqc_plots = channel.empty()
        .mix(RNA_QC_PLOTS.out.composition_pdf)
        .mix(RNA_QC_PLOTS.out.composition_jpg)
        .mix(RNA_QC_PLOTS.out.consensus_scores)
        .mix(RNA_QC_PLOTS.out.summary_verdict)
        .mix(RNA_QC_PLOTS.out.summary_verdict_group)
        .collect()

    per_biosample_status = RNA_QC_PLOTS.out.per_biosample_status
    multiqc_report = MULTIQC.out.report
}

// ============================================================================
// OUTPUT CONFIGURATION
// ============================================================================
output {
    bam_files {
        path "bam/${params.workspace}/rna/tool=star"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/bam.csv"
            header true
        }
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "rna",
             artifact: "bam",
             tool: "star",
             reference: params.genome
    }

    pca_heatmap_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/qc_plots/pca_heatmap"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "rna",
             artifact: "qc_plots_pca_heatmap"
    }

    rnaqc_summary {
        path "tables"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "rna",
             artifact: "rnaqc_summary",
             reference: params.genome
    }

    gene_count_summary {
        path "tables"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "rna",
             artifact: "gene_count_summary",
             reference: params.genome
    }

    rnaqc_metrics_combined {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/rnaqc_metrics"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "rna",
             artifact: "rnaqc_metrics_combined"
    }

    // Full-schema TSV mirroring the parquet schema (for nf-test validation)
    rnaqc_summary_tsv {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/rnaqc_metrics"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "rna",
             artifact: "rnaqc_summary_tsv"
    }

    rnaqc_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/qc_plots"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "rna",
             artifact: "rnaqc_plots"
    }

    per_biosample_status {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/index"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             artifact: "per_biosample_status"
    }

    multiqc_report {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/reports"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             artifact: "multiqc_report"
    }
}
