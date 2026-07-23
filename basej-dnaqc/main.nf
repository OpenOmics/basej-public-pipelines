nextflow.enable.dsl=2

// ============================================================================
// BASEJ-DNAQC: Single-Cell DNA QC Pipeline
// ============================================================================
// Description: Subsample to 2M reads, QC metrics, CNV cluster QC, and MultiQC report
// Outputs: Per-biosample Parquet files for Athena queries + QC report
// ============================================================================

// Import Ginkgo CNV workflow (uses proper container separation from conf/modules.config)
include { GINKO_NOPUBLISH } from './modules.nf'

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
    tuple val(sample_name), env(TOTAL_READS), env(FINAL_READS), emit: read_counts
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
        # Count both mates so FINAL_READS matches the TOTAL_READS unit (R1+R2)
        R1_SUB=\$(zcat '${sample_name}_subsampled_R1.fastq.gz' | wc -l | awk '{print int(\$1/4)}')
        if [ "${paired}" = "true" ]; then
            R2_SUB=\$(zcat '${sample_name}_subsampled_R2.fastq.gz' | wc -l | awk '{print int(\$1/4)}')
            export FINAL_READS=\$((R1_SUB + R2_SUB))
        else
            export FINAL_READS=\$R1_SUB
        fi
    fi

    echo "Final reads: \$FINAL_READS"

    # Write read counts to file for QC_PLOTS
    echo "\$TOTAL_READS" > '${sample_name}_read_counts.txt'
    echo "\$FINAL_READS" >> '${sample_name}_read_counts.txt'
    """
}

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
    tuple val(sample_name), env(TOTAL_READS), env(FINAL_READS), emit: read_counts
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
    
    # Write read counts to file for QC_PLOTS (matching SEQKIT_SAMPLE format)
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
// PROCESS: SENTIEON_ALIGN_DEDUP_METRICS
// Description: Alignment (FASTQ), deduplication, and metrics collection
//              Supports both FASTQ (BWA-MEM) and CRAM (direct processing) inputs
// ============================================================================
process SENTIEON_ALIGN_DEDUP_METRICS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(read_type)
    path fasta_ref
    path dbsnp
    path dbsnp_index
    path mills
    path mills_index
    path onekg
    path onekg_index
    val(platform)
    path base_metrics_intervals

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("*sentieonmetrics*"), emit: metrics
    tuple val(sample_name), path("${sample_name}.dedup_sentieonmetrics.txt"), emit: dedup_metrics
    path("*sentieonmetrics*"), emit: metrics_flat

    script:
    def r1 = reads[0]
    def r2 = reads.size() == 2 ? reads[1] : ""
    
    if (read_type == "CRAM") {
        // CRAM input path: Dedup + Metrics directly (no alignment)
        """
        set +u
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

        export bwt_max_mem=\$([ ${task.memory.toGiga()} -gt 30 ] && echo "30G" || echo "${task.memory.toGiga()}G")

        # LocusCollector for dedup (read from CRAM directly)
        sentieon driver -t ${task.cpus} -r '${fasta_ref}/genome.fa' -i ${r1} \\
            --algo LocusCollector --fun score_info '${sample_name}.locuscollector_score.gz'

        # Dedup with metrics (read from CRAM directly)
        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${r1} \\
            --algo Dedup --rmdup --score_info ${sample_name}.locuscollector_score.gz \\
            --metrics ${sample_name}.dedup_sentieonmetrics.txt \\
            ${sample_name}.bam

        # Collect QC metrics (read from CRAM directly)
        # QualityYield emits Q20_BASES / Q30_BASES / TOTAL_BASES — used as a Q30
        # source for CRAM input where fastp is not run (e.g. Ultima).
        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${r1} \\
            --interval ${base_metrics_intervals} \\
            --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \\
            --algo AlignmentStat ${sample_name}.alignmentstat_sentieonmetrics.txt \\
            --algo InsertSizeMetricAlgo ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \\
            --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \\
            --algo QualityYield ${sample_name}.qualityyield_sentieonmetrics.txt \\
            --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics --omit_base_output --omit_locus_stat --omit_sample_stat

        rm -f ${sample_name}.locuscollector_score.gz
        """
    } else {
        // FASTQ input path: Alignment + Dedup + Metrics (original logic)
        """
        set +u
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

        export bwt_max_mem=\$([ ${task.memory.toGiga()} -gt 30 ] && echo "30G" || echo "${task.memory.toGiga()}G")

        # BWA-MEM alignment
        sentieon bwa mem -M -Y -K 2500000000 \\
            -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:${platform}" \\
            -t ${task.cpus} '${fasta_ref}/genome.fa' '${r1}' '${r2}' | \\
            sentieon util sort -r '${fasta_ref}/genome.fa' -o '${sample_name}_sorted.bam' -t ${task.cpus} --sam2bam -i -
        
        # LocusCollector for dedup
        sentieon driver -t ${task.cpus} -r '${fasta_ref}/genome.fa' -i ${sample_name}_sorted.bam \\
            --algo LocusCollector --fun score_info '${sample_name}.locuscollector_score.gz'

        # Dedup with metrics
        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${sample_name}_sorted.bam \\
            --algo Dedup --rmdup --score_info ${sample_name}.locuscollector_score.gz \\
            --metrics ${sample_name}.dedup_sentieonmetrics.txt \\
            ${sample_name}.bam

        # Collect QC metrics
        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${sample_name}.bam \\
            --interval ${base_metrics_intervals} \\
            --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \\
            --algo AlignmentStat ${sample_name}.alignmentstat_sentieonmetrics.txt \\
            --algo InsertSizeMetricAlgo ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \\
            --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \\
            --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics --omit_base_output --omit_locus_stat --omit_sample_stat

        # Cleanup
        rm -f ${sample_name}_sorted.bam ${sample_name}.locuscollector_score.gz
        """
    }
}

// ============================================================================
// PROCESS: BWAMEM2_ALIGN
// Description: Open-source alignment using BWA-MEM2 + samtools sort
//              Equivalent to Sentieon's `bwa mem` + `util sort`
// ============================================================================
process BWAMEM2_ALIGN {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(read_type)
    path bwamem2_index_dir
    val(platform)

    output:
    tuple val(sample_name), path("${sample_name}_sorted.bam"), path("${sample_name}_sorted.bam.bai"), val(read_type), emit: bam

    script:
    def r1 = reads[0]
    def r2 = reads.size() == 2 ? reads[1] : ""

    if (read_type == "CRAM") {
        // CRAM input: decode + sort + index, no realignment (matches Sentieon CRAM path).
        // CRAM decoding requires the reference FASTA, so pipe through `samtools view
        // --reference` (using genome.fa staged in the index dir) before sorting.
        """
        set +u

        REF_FASTA="${bwamem2_index_dir}/genome.fa"
        samtools view -h -u --reference \$REF_FASTA ${r1} | \\
            samtools sort -@ ${task.cpus} -o ${sample_name}_sorted.bam -
        samtools index -@ ${task.cpus} ${sample_name}_sorted.bam
        """
    } else {
        // FASTQ input: BWA-MEM2 alignment + samtools sort.
        // Uses the pre-built BWA-MEM2 index staged in `bwamem2_index_dir`
        // (genome.fa + genome.fa.0123/.bwt.2bit.64/.amb/.ann/.pac).
        """
        set +u

        BWA_INDEX="${bwamem2_index_dir}/genome.fa"
        # BWA-MEM2 alignment piped to samtools sort
        bwa-mem2.avx2 mem -M -Y -K 2500000000 \\
            -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:${platform}" \\
            -t ${task.cpus} \$BWA_INDEX ${r1} ${r2} | \\
            samtools sort -@ ${task.cpus} -o ${sample_name}_sorted.bam -

        samtools index -@ ${task.cpus} ${sample_name}_sorted.bam
        """
    }
}

// ============================================================================
// PROCESS: SAMTOOLS_MARKDUP
// Description: Mark duplicates using samtools markdup with stats output
//              Equivalent to Sentieon's LocusCollector + Dedup
// ============================================================================
process SAMTOOLS_MARKDUP {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(sorted_bam), path(sorted_bai), val(read_type)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("${sample_name}.dedup_sentieonmetrics.txt"), emit: dedup_metrics
    path("${sample_name}.dedup_sentieonmetrics.txt"), emit: dedup_metrics_flat

    script:
    """
    set +u

    # samtools markdup requires fixmate + sort by coordinate.
    # Uncompressed (-u) intermediates skip pointless compress/decompress.
    samtools sort -n -u -@ ${task.cpus} -o ${sample_name}_namesorted.bam ${sorted_bam}
    samtools fixmate -m -u -@ ${task.cpus} ${sample_name}_namesorted.bam ${sample_name}_fixmate.bam
    samtools sort -u -@ ${task.cpus} -o ${sample_name}_resorted.bam ${sample_name}_fixmate.bam

    # Remove duplicates (-r) to match Sentieon's `Dedup --rmdup` behavior, so the
    # output BAM (and therefore preseq, coverage, and downstream metrics) is
    # consistent between the open-source and Sentieon paths. `-s` still reports the
    # duplicate counts (computed before removal) used to build the dedup metrics.
    samtools markdup -r -s -f ${sample_name}.markdup_stats.txt -@ ${task.cpus} \\
        ${sample_name}_resorted.bam ${sample_name}.bam

    samtools index -@ ${task.cpus} ${sample_name}.bam

    # Convert samtools markdup stats to Picard-style DuplicationMetrics format so
    # the QC_PLOTS Python parser (which reads READ_PAIRS_EXAMINED / READ_PAIR_DUPLICATES
    # / READ_PAIR_OPTICAL_DUPLICATES / PERCENT_DUPLICATION) works unchanged.
    #
    # samtools markdup -s emits lines like:
    #   PAIRED: <n reads>
    #   EXAMINED: <n reads>
    #   DUPLICATE PAIR: <n reads>
    #   DUPLICATE PAIR OPTICAL: <n reads>
    # Picard counts PAIRS, so divide read counts by 2.
    get_stat() { grep -E "^\$1:" ${sample_name}.markdup_stats.txt | head -1 | awk -F': ' '{print \$2+0}'; }

    PAIRED_READS=\$(get_stat "PAIRED")
    DUP_PAIR_READS=\$(get_stat "DUPLICATE PAIR")
    DUP_PAIR_OPT_READS=\$(get_stat "DUPLICATE PAIR OPTICAL")
    UNPAIRED_READS=\$(get_stat "SINGLE")
    UNPAIRED_DUPS=\$(get_stat "DUPLICATE SINGLE")
    UNMAPPED=\$(get_stat "EXCLUDED")

    # Convert read counts to pair counts (Picard convention)
    READ_PAIRS_EXAMINED=\$(awk -v r=\$PAIRED_READS 'BEGIN{printf "%d", r/2}')
    READ_PAIR_DUPLICATES=\$(awk -v r=\$DUP_PAIR_READS 'BEGIN{printf "%d", r/2}')
    READ_PAIR_OPTICAL_DUPLICATES=\$(awk -v r=\$DUP_PAIR_OPT_READS 'BEGIN{printf "%d", r/2}')

    # PERCENT_DUPLICATION = (pair dups + unpaired dups) / (pairs examined + unpaired examined)
    PCT_DUP=\$(awk -v dp=\$READ_PAIR_DUPLICATES -v du=\$UNPAIRED_DUPS -v ep=\$READ_PAIRS_EXAMINED -v eu=\$UNPAIRED_READS \\
        'BEGIN { denom = ep + eu; if (denom > 0) printf "%.6f", (dp + du) / denom; else printf "0" }')

    cat > ${sample_name}.dedup_sentieonmetrics.txt << EOF
## METRICS CLASS	picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
${sample_name}	\${UNPAIRED_READS}	\${READ_PAIRS_EXAMINED}	0	\${UNMAPPED}	\${UNPAIRED_DUPS}	\${READ_PAIR_DUPLICATES}	\${READ_PAIR_OPTICAL_DUPLICATES}	\${PCT_DUP}	0
EOF

    # Cleanup intermediate files
    rm -f ${sample_name}_namesorted.bam ${sample_name}_fixmate.bam ${sample_name}_resorted.bam
    """
}

// ============================================================================
// PROCESS: PICARD_METRICS
// Description: Collect alignment/insert/GC/quality metrics with Picard tools
//              Plus per-interval coverage via samtools depth + awk to mimic
//              Sentieon's CoverageMetrics sample_interval_summary format.
//              Equivalent to Sentieon's GCBias + AlignmentStat + InsertSizeMetricAlgo
//                            + QualityYield + MeanQualityByCycle + CoverageMetrics
// ============================================================================
process PICARD_METRICS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path fasta_ref
    path base_metrics_intervals

    output:
    tuple val(sample_name), path("*sentieonmetrics*"), emit: metrics
    path("*sentieonmetrics*"), emit: metrics_flat

    script:
    """
    set +u

    # Picard metric suite — single pass via CollectMultipleMetrics. Runs
    # AlignmentSummary, InsertSize, GcBias, MeanQualityByCycle and QualityYield
    # as PROGRAM modules in ONE pass over the dedup BAM instead of five separate
    # full-BAM passes. Same Picard collectors -> identical values; only the file
    # names differ, so we rename to the *sentieonmetrics* names the parser and
    # MultiQC expect (MultiQC detects by the internal "## METRICS CLASS picard"
    # header, not the filename).
    picard CollectMultipleMetrics \\
        R=${fasta_ref}/genome.fa \\
        I=${bam} \\
        O=${sample_name}.multimetrics \\
        PROGRAM=null \\
        PROGRAM=CollectAlignmentSummaryMetrics \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=MeanQualityByCycle \\
        PROGRAM=CollectQualityYieldMetrics

    # Rename single-pass outputs to the expected *sentieonmetrics* names. Some
    # modules emit no file when there is no eligible data (e.g. InsertSize with
    # zero mapped pairs), so fall back to an empty placeholder the parser skips.
    mv_or_touch() { if [ -f "\$1" ]; then mv "\$1" "\$2"; else touch "\$2"; fi; }
    mv_or_touch ${sample_name}.multimetrics.insert_size_metrics      ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.gc_bias.summary_metrics  ${sample_name}.gcbias_summary.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.gc_bias.detail_metrics   ${sample_name}.gcbias.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.quality_by_cycle_metrics ${sample_name}.meanqualitybycycle.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.quality_yield_metrics    ${sample_name}.qualityyield_sentieonmetrics.txt

    # AlignmentSummaryMetrics: drop Picard's 3 trailing columns (SAMPLE, LIBRARY,
    # READ_GROUP) so the CATEGORY data row and header have matching field counts
    # in the parser (otherwise align_* / PCT_CHIMERAS come out empty).
    # Tolerates a missing file (empty input).
    python3 << 'PYEOF'
import os
inp = "${sample_name}.multimetrics.alignment_summary_metrics"
out = "${sample_name}.alignmentstat_sentieonmetrics.txt"
cats = ("CATEGORY", "FIRST_OF_PAIR", "SECOND_OF_PAIR", "PAIR", "UNPAIRED")
if not os.path.exists(inp):
    open(out, "w").close()
else:
    with open(inp) as fh, open(out, "w") as of:
        for line in fh:
            if line.startswith("#") or not line.strip():
                of.write(line)
                continue
            if line.split("\\t")[0] in cats:
                parts = line.rstrip("\\n").split("\\t")
                # drop the last 3 columns (SAMPLE, LIBRARY, READ_GROUP)
                of.write("\\t".join(parts[:-3]) + "\\n")
            else:
                of.write(line)
PYEOF
    rm -f ${sample_name}.multimetrics.alignment_summary_metrics

    # Per-interval coverage via samtools depth + awk → cov_*.sample_interval_summary.
    # The parser only reads column 2 (per-chromosome summed depth) to compute
    # cov_total_bases / cov_chrm_bases / pct_chrm, so `samtools depth` WITHOUT
    # `-a` gives identical sums (zero-depth positions contribute nothing) while
    # skipping the costly all-positions output.
    samtools depth -b ${base_metrics_intervals} ${bam} > ${sample_name}.depth.tsv

    awk 'BEGIN {
        OFS="\\t"
        print "Target","total_coverage","average_coverage","sample_name_total_cvg","sample_name_mean_cvg","sample_name_granular_Q1","sample_name_median","sample_name_granular_Q3","sample_name_pct_above_15"
    }
    {
        key=\$1
        chr_total[key] += \$3
        chr_count[key] += 1
        if (\$3 >= 15) chr_above15[key] += 1
    }
    END {
        for (k in chr_total) {
            mean = (chr_count[k] > 0) ? chr_total[k]/chr_count[k] : 0
            pct15 = (chr_count[k] > 0) ? (chr_above15[k]/chr_count[k])*100 : 0
            print k, chr_total[k], mean, chr_total[k], mean, 0, 0, 0, pct15
        }
    }' ${sample_name}.depth.tsv > ${sample_name}.cov_sentieonmetrics.sample_interval_summary

    rm -f ${sample_name}.depth.tsv ${sample_name}.multimetrics.*.pdf
    """
}

// ============================================================================
// PROCESS: PRESEQ
// Description: Library complexity estimation
// ============================================================================
process PRESEQ {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    tuple val(sample_name), path("${sample_name}_preseq.txt"), emit: complexity
    path("${sample_name}_preseq.txt"), emit: complexity_flat

    script:
    """
    set +u

    echo "0" > ${sample_name}_preseq.txt
    {
        bam2mr -seg_len 100000 -o ${sample_name}.mr ${bam} 2>/dev/null && \\
        preseq gc_extrap -o ${sample_name}.curve ${sample_name}.mr && \\
        tail -n 1 ${sample_name}.curve | cut -f 2 | awk '{printf "%.0f\\n", \$1}' > ${sample_name}_preseq.txt
    } || {
        echo "Preseq failed (likely low complexity library), using 0"
    }
    
    echo "Estimated library complexity: \$(cat ${sample_name}_preseq.txt)"
    """
}

// ============================================================================
// ============================================================================
// PROCESS: QC_PLOTS (CONSOLIDATED)
// Description: Unified process combining:
//   1. DNAQC_METRICS_TO_PARQUET - Parse per-sample metrics, create Parquets + TSVs
//   2. MERGE_QC_METRICS - Merge all TSV files into combined_selected_metrics.txt
//   3. QC_PLOTS - Generate visualizations and consensus scores
// 
// This consolidated process:
//   - Accepts collected metric tuples from all samples (via glob in script)
//   - Parses all metrics (Sentieon, Fastp, Preseq, Ginkgo) per sample
//   - Writes individual Parquet files for each sample
//   - Merges all TSV metrics into single combined file
//   - Generates CNV summary, plots, and consensus scores
// ============================================================================
process QC_PLOTS {   
    tag "qc_plots"

    input:
    path(rds_file)
    path(seg_copy_file)
    path(ginkgo_metrics)
    path(metadata_file)
    path(plot_qc_config)
    path(sentieon_metrics)
    path(fastp_jsons)
    path(preseq_files)
    path(read_counts_files)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(user)

    output:
    path("dnaqc_summary/workspace=*/workflow_id=*/biosample=*/output.parquet"), emit: parquet
    path("*_selected_metrics_mqc.txt"), emit: mqc_metrics
    path("AllSample-GinkgoSegmentSummary.txt")
    path("combined_selected_metrics.txt"), emit: combined_metrics
    path("dnaqc_all_metrics.tsv"), emit: summary_tsv
    path("nf-preseq-pipeline_all_metrics_mqc.txt"), emit: allmetrics_with_cnv
    path("DNA-QC_ConsensusScores.txt"), emit: consensus_scores
    path("DNA-QC_ConsensusScores_SummaryTable_mqc.txt"), emit: consensus_summary
    path("ConsensusScores_SummaryTableByGroup.txt"), emit: consensus_group_summary
    path("QC_composition.pdf"), emit: composition_pdf
    path("CNV-Quadrants.pdf"), emit: cnv_quadrants_pdf
    path("QC_composition_mqc.jpg"), emit: composition_jpg
    path("CNV-Quadrants_mqc.jpg"), emit: cnv_quadrants_jpg
    path("per_biosample_status.csv"), emit: per_biosample_status

    script:
    def has_ginkgo = ginkgo_metrics.name != 'NO_GINKGO_METRICS'
    """
    #!/bin/bash
    set -euo pipefail

    echo "=============================================="
    echo "QC_PLOTS: Starting consolidated QC workflow"
    echo "  1. Parse per-sample metrics → Parquets + TSVs"
    echo "  2. Merge all metrics into combined file"
    echo "  3. Generate visualizations"
    echo "=============================================="
    echo ""

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
has_ginkgo = ${has_ginkgo ? 'True' : 'False'}

# Discover all unique samples from Sentieon metric files
sample_names = set()
for f in glob.glob("*.alignmentstat_sentieonmetrics.txt"):
    # Extract sample name: sampleA.alignmentstat_sentieonmetrics.txt -> sampleA
    basename = os.path.basename(f)
    sample = basename.split(".alignmentstat")[0]
    sample_names.add(sample)

if not sample_names:
    print("ERROR: No sample metric files found")
    sys.exit(1)

print(f"Found {len(sample_names)} samples: {sorted(sample_names)}")

# Load Ginkgo metrics if available
df_ginkgo = None
if has_ginkgo:
    try:
        df_ginkgo = pd.read_csv("${ginkgo_metrics}", sep="\\t")
        if "SampleId" in df_ginkgo.columns:
            df_ginkgo = df_ginkgo.rename(columns={"SampleId": "biosample"})
        df_ginkgo["biosample"] = df_ginkgo["biosample"].str.replace("_sorted", "", regex=False)
        print(f"Loaded Ginkgo metrics for {len(df_ginkgo)} samples")
    except Exception as e:
        print(f"Warning: Could not load Ginkgo metrics: {e}")
        df_ginkgo = None

# Load Ginkgo SegCopy if available
df_segcopy = None
if has_ginkgo:
    try:
        df_segcopy = pd.read_csv("${seg_copy_file}", sep="\\t")
        print(f"Loaded SegCopy data")
    except Exception as e:
        print(f"Warning: Could not load SegCopy: {e}")

all_summaries = []

# Process each sample
for sample in sorted(sample_names):
    summary = {
        "biosample": sample,
        "dataset_id": dataset_id,
        "pipeline": "bj-dnaqc",
        "pipeline_version": pipeline_version,
        "molecule_type": "dna",
        "user": user,
    }
    
    # Find per-sample metric files (use glob with sample name)
    fastp_files = glob.glob(f"{sample}*fastp.json")
    preseq_files = glob.glob(f"{sample}*preseq.txt")
    read_count_files = glob.glob(f"{sample}*read_counts.txt")
    
    # ========== Parse Read Counts ==========
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
        except Exception as e:
            print(f"Warning: Could not parse fastp for {sample}: {e}")
    
    # ========== Parse Alignment Stats ==========
    # Picard CollectAlignmentSummaryMetrics / Sentieon AlignmentStat emits one
    # row per CATEGORY: FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR (paired-end), or
    # UNPAIRED (single-end, e.g. Ultima CRAM). Prefer PAIR; fall back to
    # UNPAIRED so single-end data still populates total_reads, mismatch rate,
    # mean_read_length, etc. PCT_CHIMERAS is undefined for single-end and
    # remains 0.0 in that case.
    for f in glob.glob(f"{sample}*.alignmentstat_sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                headers = []
                rows_by_cat = {}
                for line in lines:
                    if line.startswith("CATEGORY"):
                        headers = line.strip().split("\\t")
                    elif headers and line.strip() and not line.startswith("#"):
                        parts = line.strip().split("\\t")
                        if len(parts) == len(headers):
                            cat = parts[0]
                            if cat in ("PAIR", "UNPAIRED"):
                                rows_by_cat[cat] = dict(zip(headers, parts))
                d = rows_by_cat.get("PAIR") or rows_by_cat.get("UNPAIRED")
                if d:
                    summary["align_total_reads"] = int(d.get("TOTAL_READS", 0))
                    summary["align_pf_reads"] = int(d.get("PF_READS", 0))
                    summary["align_pf_reads_aligned"] = int(d.get("PF_READS_ALIGNED", 0))
                    summary["align_pf_hq_aligned_reads"] = int(d.get("PF_HQ_ALIGNED_READS", 0))
                    summary["align_pf_hq_aligned_bases"] = int(d.get("PF_HQ_ALIGNED_BASES", 0))
                    summary["align_pf_hq_aligned_q20_bases"] = int(d.get("PF_HQ_ALIGNED_Q20_BASES", 0))
                    summary["align_pf_mismatch_rate"] = float(d.get("PF_MISMATCH_RATE", 0))
                    summary["align_pf_hq_error_rate"] = float(d.get("PF_HQ_ERROR_RATE", 0))
                    summary["align_pf_indel_rate"] = float(d.get("PF_INDEL_RATE", 0))
                    summary["align_mean_read_length"] = float(d.get("MEAN_READ_LENGTH", 0))
                    summary["align_reads_aligned_in_pairs"] = int(d.get("READS_ALIGNED_IN_PAIRS", 0))
                    summary["align_pct_reads_aligned_in_pairs"] = float(d.get("PCT_READS_ALIGNED_IN_PAIRS", 0))
                    summary["align_bad_cycles"] = int(d.get("BAD_CYCLES", 0))
                    summary["align_strand_balance"] = float(d.get("STRAND_BALANCE", 0))
                    summary["align_pct_chimeras"] = float(d.get("PCT_CHIMERAS", 0)) if d.get("PCT_CHIMERAS") else 0
                    summary["align_pct_adapter"] = float(d.get("PCT_ADAPTER", 0)) if d.get("PCT_ADAPTER") else 0
                    total = summary["align_total_reads"]
                    aligned = summary["align_pf_reads_aligned"]
                    summary["pct_aligned"] = (aligned / total) if total > 0 else 0
                    summary["pct_pf"] = (summary["align_pf_reads"] / total) if total > 0 else 0
                    summary["pct_error"] = summary["align_pf_mismatch_rate"]
                    summary["pct_chimeras"] = summary["align_pct_chimeras"]
        except Exception as e:
            print(f"Warning: Could not parse alignment stats for {sample}: {e}")

    # ========== Parse Sentieon QualityYield (CRAM input only) ==========
    # When fastp is not run (CRAM path, e.g. Ultima), derive Q30 rate from
    # Sentieon QualityYield. Falls back gracefully if the file is absent.
    for f in glob.glob(f"{sample}*.qualityyield_sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                headers = []
                for i, line in enumerate(lines):
                    if line.startswith("TOTAL_READS") or line.startswith("TOTAL_BASES"):
                        headers = line.strip().split("\\t")
                        if i + 1 < len(lines):
                            parts = lines[i + 1].strip().split("\\t")
                            if len(parts) == len(headers):
                                d = dict(zip(headers, parts))
                                total_bases = float(d.get("TOTAL_BASES", 0) or 0)
                                q20_bases = float(d.get("Q20_BASES", 0) or 0)
                                q30_bases = float(d.get("Q30_BASES", 0) or 0)
                                summary["sentieon_total_bases"] = total_bases
                                summary["sentieon_q20_rate"] = (q20_bases / total_bases) if total_bases > 0 else 0
                                summary["sentieon_q30_rate"] = (q30_bases / total_bases) if total_bases > 0 else 0
                                # Populate fastp_q30_rate from Sentieon when fastp
                                # was not run (CRAM path); preserves the existing
                                # MultiQC column / parquet schema.
                                if not summary.get("fastp_q30_rate"):
                                    summary["fastp_q30_rate"] = summary["sentieon_q30_rate"]
                        break
        except Exception as e:
            print(f"Warning: Could not parse QualityYield for {sample}: {e}")

    
    # ========== Parse Dedup Metrics ==========
    for f in glob.glob(f"{sample}*.dedup_sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                headers = []
                for i, line in enumerate(lines):
                    if line.startswith("LIBRARY"):
                        headers = line.strip().split("\\t")
                        if i+1 < len(lines):
                            parts = lines[i+1].strip().split("\\t")
                            d = dict(zip(headers, parts))
                            summary["dedup_unpaired_reads_examined"] = int(d.get("UNPAIRED_READS_EXAMINED", 0))
                            summary["dedup_read_pairs_examined"] = int(d.get("READ_PAIRS_EXAMINED", 0))
                            summary["dedup_secondary_or_supplementary_reads"] = int(d.get("SECONDARY_OR_SUPPLEMENTARY_RDS", 0))
                            summary["dedup_unmapped_reads"] = int(d.get("UNMAPPED_READS", 0))
                            summary["dedup_unpaired_read_duplicates"] = int(d.get("UNPAIRED_READ_DUPLICATES", 0))
                            summary["dedup_read_pair_duplicates"] = int(d.get("READ_PAIR_DUPLICATES", 0))
                            summary["dedup_read_pair_optical_duplicates"] = int(d.get("READ_PAIR_OPTICAL_DUPLICATES", 0))
                            summary["dedup_pct_duplication"] = float(d.get("PERCENT_DUPLICATION", 0))
                            summary["dedup_estimated_library_size"] = int(d.get("ESTIMATED_LIBRARY_SIZE", 0)) if d.get("ESTIMATED_LIBRARY_SIZE") else 0
                            summary["pct_duplication"] = summary["dedup_pct_duplication"]
                            read_pairs_examined = summary["dedup_read_pairs_examined"]
                            optical_dups = summary["dedup_read_pair_optical_duplicates"]
                            total_dups = summary["dedup_read_pair_duplicates"]
                            summary["pct_optical_duplicates"] = (optical_dups / read_pairs_examined) if read_pairs_examined > 0 else 0
                            summary["pct_pcr_duplicates"] = ((total_dups - optical_dups) / read_pairs_examined) if read_pairs_examined > 0 else 0
                        break
        except Exception as e:
            print(f"Warning: Could not parse dedup for {sample}: {e}")
    
    # ========== Parse Insert Size Metrics ==========
    for f in glob.glob(f"{sample}*.insertsizemetricalgo.sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                headers = []
                for i, line in enumerate(lines):
                    if line.startswith("MEDIAN_INSERT_SIZE"):
                        headers = line.strip().split("\\t")
                        if i+1 < len(lines):
                            parts = lines[i+1].strip().split("\\t")
                            d = dict(zip(headers, parts))
                            summary["insert_median"] = float(d.get("MEDIAN_INSERT_SIZE", 0))
                            summary["insert_mode"] = float(d.get("MODE_INSERT_SIZE", 0)) if d.get("MODE_INSERT_SIZE") else 0
                            summary["insert_median_absolute_deviation"] = float(d.get("MEDIAN_ABSOLUTE_DEVIATION", 0))
                            summary["insert_min"] = float(d.get("MIN_INSERT_SIZE", 0)) if d.get("MIN_INSERT_SIZE") else 0
                            summary["insert_max"] = float(d.get("MAX_INSERT_SIZE", 0)) if d.get("MAX_INSERT_SIZE") else 0
                            summary["insert_mean"] = float(d.get("MEAN_INSERT_SIZE", 0))
                            summary["insert_std"] = float(d.get("STANDARD_DEVIATION", 0))
                            summary["insert_read_pairs"] = int(d.get("READ_PAIRS", 0)) if d.get("READ_PAIRS") else 0
                            summary["insert_pair_orientation"] = d.get("PAIR_ORIENTATION", "")
                            summary["insert_width_of_10_pct"] = float(d.get("WIDTH_OF_10_PERCENT", 0)) if d.get("WIDTH_OF_10_PERCENT") else 0
                            summary["insert_width_of_20_pct"] = float(d.get("WIDTH_OF_20_PERCENT", 0)) if d.get("WIDTH_OF_20_PERCENT") else 0
                            summary["insert_width_of_30_pct"] = float(d.get("WIDTH_OF_30_PERCENT", 0)) if d.get("WIDTH_OF_30_PERCENT") else 0
                            summary["insert_width_of_40_pct"] = float(d.get("WIDTH_OF_40_PERCENT", 0)) if d.get("WIDTH_OF_40_PERCENT") else 0
                            summary["insert_width_of_50_pct"] = float(d.get("WIDTH_OF_50_PERCENT", 0)) if d.get("WIDTH_OF_50_PERCENT") else 0
                            summary["insert_width_of_60_pct"] = float(d.get("WIDTH_OF_60_PERCENT", 0)) if d.get("WIDTH_OF_60_PERCENT") else 0
                            summary["insert_width_of_70_pct"] = float(d.get("WIDTH_OF_70_PERCENT", 0)) if d.get("WIDTH_OF_70_PERCENT") else 0
                            summary["insert_width_of_80_pct"] = float(d.get("WIDTH_OF_80_PERCENT", 0)) if d.get("WIDTH_OF_80_PERCENT") else 0
                            summary["insert_width_of_90_pct"] = float(d.get("WIDTH_OF_90_PERCENT", 0)) if d.get("WIDTH_OF_90_PERCENT") else 0
                            summary["insert_width_of_95_pct"] = float(d.get("WIDTH_OF_95_PERCENT", 0)) if d.get("WIDTH_OF_95_PERCENT") else 0
                            summary["insert_width_of_99_pct"] = float(d.get("WIDTH_OF_99_PERCENT", 0)) if d.get("WIDTH_OF_99_PERCENT") else 0
                            summary["insert_size"] = summary["insert_median"]
                        break
        except Exception as e:
            print(f"Warning: Could not parse insert size for {sample}: {e}")
    
    # ========== Parse GC Bias Summary ==========
    for f in glob.glob(f"{sample}*.gcbias_summary.sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                headers = []
                for i, line in enumerate(lines):
                    if line.startswith("ACCUMULATION"):
                        headers = line.strip().split("\\t")
                        if i+1 < len(lines):
                            parts = lines[i+1].strip().split("\\t")
                            d = dict(zip(headers, parts))
                            summary["gc_total_clusters"] = int(d.get("TOTAL_CLUSTERS", 0)) if d.get("TOTAL_CLUSTERS") else 0
                            summary["gc_aligned_reads"] = int(d.get("ALIGNED_READS", 0)) if d.get("ALIGNED_READS") else 0
                            summary["gc_at_dropout"] = float(d.get("AT_DROPOUT", 0))
                            summary["gc_gc_dropout"] = float(d.get("GC_DROPOUT", 0))
                            summary["gc_gc_nc_0_19"] = float(d.get("GC_NC_0_19", 0)) if d.get("GC_NC_0_19") else 0
                            summary["gc_gc_nc_20_39"] = float(d.get("GC_NC_20_39", 0)) if d.get("GC_NC_20_39") else 0
                            summary["gc_gc_nc_40_59"] = float(d.get("GC_NC_40_59", 0)) if d.get("GC_NC_40_59") else 0
                            summary["gc_gc_nc_60_79"] = float(d.get("GC_NC_60_79", 0)) if d.get("GC_NC_60_79") else 0
                            summary["gc_gc_nc_80_100"] = float(d.get("GC_NC_80_100", 0)) if d.get("GC_NC_80_100") else 0
                        break
        except Exception as e:
            print(f"Warning: Could not parse GC bias for {sample}: {e}")
    
    # ========== Parse Coverage Metrics ==========
    for f in glob.glob(f"{sample}*.cov_sentieonmetrics.sample_interval_summary"):
        try:
            total_bases, chrm_bases = 0, 0
            with open(f, 'r') as fh:
                for line in fh:
                    if line.startswith("#") or line.startswith("Target"):
                        continue
                    parts = line.strip().split("\\t")
                    if len(parts) >= 2:
                        chrom = parts[0].split(":")[0]
                        try:
                            bases = int(float(parts[1]))
                        except (ValueError, IndexError):
                            bases = 0
                        total_bases += bases
                        if chrom in ["chrM", "MT"]:
                            chrm_bases = bases
            summary["cov_total_bases"] = total_bases
            summary["cov_chrm_bases"] = chrm_bases
            summary["pct_chrm"] = (chrm_bases / total_bases) if total_bases > 0 else 0
        except Exception as e:
            print(f"Warning: Could not parse coverage for {sample}: {e}")
    
    # ========== Parse Preseq (library complexity) ==========
    if preseq_files:
        try:
            with open(preseq_files[0], 'r') as f:
                summary["preseq_count"] = int(f.read().strip())
        except:
            summary["preseq_count"] = 0
    else:
        summary["preseq_count"] = 0
    
    # ========== Initialize Ginkgo fields ==========
    summary["ginkgo_cnv_mapd"] = None
    summary["ginkgo_cnv_skew"] = None
    summary["ginkgo_ploidy"] = None
    summary["ginkgo_average_ploidy"] = None
    
    # ========== Parse Ginkgo metrics (ploidy only — MAPD/SKEW come from cnvSummarizer in Step 4) ==========
    if has_ginkgo and df_ginkgo is not None:
        try:
            sample_row = df_ginkgo[df_ginkgo["biosample"] == sample]
            if len(sample_row) > 0:
                row = sample_row.iloc[0]
                summary["ginkgo_ploidy"] = float(row.get("GenomePloidy")) if pd.notna(row.get("GenomePloidy")) else None
        except Exception as e:
            print(f"Warning: Could not parse Ginkgo summary for {sample}: {e}")
    
    # ========== Parse SegCopy for average ploidy ==========
    if has_ginkgo and df_segcopy is not None:
        try:
            sample_col = None
            for col in df_segcopy.columns:
                if col.replace("_sorted", "") == sample:
                    sample_col = col
                    break
            if sample_col:
                summary["ginkgo_average_ploidy"] = float(df_segcopy[sample_col].mean())
        except Exception as e:
            print(f"Warning: Could not parse SegCopy for {sample}: {e}")
    
    all_summaries.append(summary)

    # ========== Save summary as JSON (parquets written after R scoring in Step 8) ==========
    with open(f"{sample}_summary.json", 'w') as jf:
        json.dump(summary, jf, default=str)

    # ========== Write individual TSV for MultiQC/R input (QC_Status filled in Step 8) ==========
    subsampled = "Yes" if summary.get('subsampled') else "No"

    with open(f"{sample}_selected_metrics_mqc.txt", 'w') as f:
        f.write("# id: 'dnaqc_summary'\\n")
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'QC Summary'\\n")
        f.write("# description: 'Per-sample alignment, CNV noise, and library complexity QC metrics with pass/fail status.'\\n")
        f.write("# pconfig:\\n")
        f.write("#   id: 'dnaqc_summary_table'\\n")
        f.write("sample_name\\tpreseq_count\\tPCT_CHIMERAS\\tchrM\\tTotalReads\\tFinalReads\\tinsert_size\\tfastp_q30_rate\\n")
        f.write(f"{sample}\\t")
        f.write(f"{summary.get('preseq_count', 0)}\\t{summary.get('align_pct_chimeras', 0):.4f}\\t")
        f.write(f"{summary.get('pct_chrm', 0):.4f}\\t")
        f.write(f"{total_reads:,}\\t{final_reads:,}\\t")
        f.write(f"{summary.get('insert_size', 0):.0f}\\t{summary.get('fastp_q30_rate', 0):.4f}\\n")

    print(f"  {sample}: metrics parsed, QC status pending R scoring")

print(f"Processed {len(all_summaries)} samples: created JSONs and TSVs")
PYEOF

    # ========== SECTION 2: Merge all per-sample TSV files ==========
    echo "Step 2: Merging per-sample metrics into combined file..."
    
    python3 << 'PYEOF'
import glob

tsv_files = sorted(glob.glob("*_selected_metrics_mqc.txt"))
if not tsv_files:
    print("ERROR: No metrics files found")
    exit(1)

with open("combined_selected_metrics.txt", 'w') as outf:
    header_written = False
    for tsv_file in tsv_files:
        with open(tsv_file, 'r') as inf:
            non_comment_count = 0
            for line in inf:
                if line.startswith("#"):
                    continue
                non_comment_count += 1
                if not header_written:
                    outf.write(line)
                    header_written = True
                elif non_comment_count > 1:
                    outf.write(line)

print(f"Merged {len(tsv_files)} sample metrics files")
PYEOF

    # ========== SECTION 3: Validate input files and run R scripts ==========
    echo "Step 3: Validating input files for visualization..."
    
    if [[ ! -f "${rds_file}" ]]; then
        echo "ERROR: RDS file not found: ${rds_file}"
        exit 1
    fi
    echo "  ✓ RDS file: ${rds_file}"
    
    if [[ ! -f "${seg_copy_file}" ]]; then
        echo "ERROR: SegCopy file not found: ${seg_copy_file}"
        exit 1
    fi
    echo "  ✓ SegCopy file: ${seg_copy_file}"
    
    if [[ ! -f "${metadata_file}" ]]; then
        echo "ERROR: Metadata file not found: ${metadata_file}"
        exit 1
    fi
    echo "  ✓ Metadata file: ${metadata_file}"
    
    if [[ ! -f "${plot_qc_config}" ]]; then
        echo "ERROR: Plot QC config not found: ${plot_qc_config}"
        exit 1
    fi
    echo "  ✓ Config file: ${plot_qc_config}"
    
    echo ""
    echo "Step 4: Running CNV Summarizer..."
    Rscript /usr/local/bin/cnvSummarizer.R \\
        --rds_file ${rds_file} \\
        --out_file AllSample-GinkgoSegmentSummary.txt
    echo "  ✓ Finished cnvSummarizer.R"

    echo ""
    echo "Step 5: Merging metrics with CNV summary..."
    python3 << 'PYEOF'
import pandas as pd
import numpy as np

df_metrics = pd.read_csv("combined_selected_metrics.txt", sep="\\t")
df_cnv = pd.read_csv("AllSample-GinkgoSegmentSummary.txt", sep="\\t")

df_cnv_subset = df_cnv[["SampleId", "MAPD_CNV_Log2", "SKEW_CNV"]].copy()
df_cnv_subset.rename(columns={"SampleId": "sample_name"}, inplace=True)

# Replace "NA" strings with NaN for proper handling
df_cnv_subset["SKEW_CNV"] = pd.to_numeric(df_cnv_subset["SKEW_CNV"], errors='coerce')

df_merged = df_metrics.merge(df_cnv_subset, on="sample_name", how="left")
df_merged["SKEW_CNV"] = df_merged["SKEW_CNV"].fillna(999)

df_merged.to_csv("nf-preseq-pipeline_all_metrics_mqc_withcnv.txt", sep="\\t", index=False)
print(f"Merged metrics with CNV data")
PYEOF
    echo "  ✓ Finished merging"

    echo ""
    echo "Step 6: Generating QC composition plots..."
    Rscript /usr/local/bin/function_plot_qc_dna.R \\
        --seg_copy_file ${seg_copy_file} \\
        --metrics_file nf-preseq-pipeline_all_metrics_mqc_withcnv.txt \\
        --metadata_file ${metadata_file} \\
        --cnv_summary_file AllSample-GinkgoSegmentSummary.txt \\
        --plot_qc_config ${plot_qc_config}
    echo "  ✓ Finished function_plot_qc_dna.R"

    echo ""
    echo "Step 7: Generating CNV quadrants plot..."
    Rscript /usr/local/bin/function_cnv_quadrants_qc.R \\
        --metrics_file nf-preseq-pipeline_all_metrics_mqc_withcnv.txt
    echo "  ✓ Finished function_cnv_quadrants_qc.R"
    
    # Rename the _withcnv file to final output name (matches original workflow)
    mv nf-preseq-pipeline_all_metrics_mqc_withcnv.txt nf-preseq-pipeline_all_metrics_mqc.txt

    echo ""
    echo "Step 8: Writing parquets and finalising QC status from R consensus scores..."
    python3 << 'PYEOF'
import os, json, glob
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

workspace  = "${workspace}"
workflow_id = "${workflow_id}"

# ========== Load MAPD/SKEW from cnvSummarizer output (written in Step 4) ==========
# AllSample-GinkgoSegmentSummary.txt is overwritten by cnvSummarizer.R in Step 4
# and contains MAPD_CNV_Log2 and SKEW_CNV columns
cnv_mapd_by_sample = {}
cnv_skew_by_sample = {}
try:
    df_cnvsumm = pd.read_csv("AllSample-GinkgoSegmentSummary.txt", sep="\\t")
    for _, row in df_cnvsumm.iterrows():
        sid = str(row.get("SampleId", "")).replace("_sorted", "")
        mapd_val = row.get("MAPD_CNV_Log2")
        skew_val = row.get("SKEW_CNV")
        cnv_mapd_by_sample[sid] = float(mapd_val) if pd.notna(mapd_val) else None
        cnv_skew_by_sample[sid] = float(skew_val) if pd.notna(skew_val) else None
except Exception as e:
    print(f"Warning: Could not load cnvSummarizer MAPD/SKEW: {e}")

# ========== Read R consensus scores ==========
df_scores = pd.read_csv("DNA-QC_ConsensusScores.txt", sep="\\t")
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

# ========== Schema column lists ==========
_str_cols    = ['biosample','dataset_id','pipeline','pipeline_version','molecule_type',
                'insert_pair_orientation','ginkgo_cnv_mapd','ginkgo_cnv_skew','ginkgo_ploidy',
                'qc_status','workflow_id','workspace','user']
_double_cols = ['fastp_q20_rate','fastp_q30_rate','fastp_gc_content','fastp_q30_rate_after',
                'fastp_duplication_rate','align_pf_mismatch_rate','align_pf_hq_error_rate',
                'align_pf_indel_rate','align_mean_read_length','align_pct_reads_aligned_in_pairs',
                'align_strand_balance','align_pct_chimeras','align_pct_adapter',
                'sentieon_q20_rate','sentieon_q30_rate',
                'pct_aligned','pct_pf','pct_error','pct_chimeras',
                'pct_optical_duplicates','pct_pcr_duplicates','dedup_pct_duplication','pct_duplication','pct_chrm',
                'insert_median','insert_median_absolute_deviation','insert_min','insert_max',
                'insert_mean','insert_std','insert_width_of_10_pct','insert_width_of_20_pct',
                'insert_width_of_30_pct','insert_width_of_40_pct','insert_width_of_50_pct',
                'insert_width_of_60_pct','insert_width_of_70_pct','insert_width_of_80_pct',
                'insert_width_of_90_pct','insert_width_of_99_pct','insert_size',
                'gc_at_dropout','gc_gc_dropout','gc_gc_nc_0_19','gc_gc_nc_20_39',
                'gc_gc_nc_40_59','gc_gc_nc_60_79','gc_gc_nc_80_100','ginkgo_average_ploidy']
_bigint_cols = ['total_reads','final_reads','fastp_total_reads','fastp_total_bases',
                'fastp_read1_mean_length','fastp_read2_mean_length','fastp_reads_after_filter',
                'fastp_adapter_trimmed_reads','fastp_adapter_trimmed_bases',
                'align_total_reads','align_pf_reads','align_pf_reads_aligned',
                'align_pf_hq_aligned_reads','align_pf_hq_aligned_bases','align_pf_hq_aligned_q20_bases',
                'align_reads_aligned_in_pairs','align_bad_cycles',
                'dedup_unpaired_reads_examined','dedup_read_pairs_examined',
                'dedup_secondary_or_supplementary_reads','dedup_unmapped_reads',
                'dedup_unpaired_read_duplicates','dedup_read_pair_duplicates',
                'dedup_read_pair_optical_duplicates','dedup_estimated_library_size',
                'insert_mode','insert_read_pairs','insert_width_of_95_pct',
                'gc_total_clusters','gc_aligned_reads','cov_total_bases','cov_chrm_bases',
                'sentieon_total_bases',
                'preseq_count','qc_score']
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

    # Override score for zero-read samples — R scoring assigns non-zero scores
    # to samples with all-zero metrics, but they should show as FAIL / NA
    if int(summary.get("total_reads", 0) or 0) == 0:
        composite_score = None
        qc_status = "FAIL"

    summary["qc_status"] = qc_status
    summary["qc_score"]  = int(composite_score) if pd.notna(composite_score) else None

    # Update MAPD/SKEW from cnvSummarizer output (Step 4) — the JSON has None from Step 1
    summary["ginkgo_cnv_mapd"] = cnv_mapd_by_sample.get(sample)
    summary["ginkgo_cnv_skew"] = cnv_skew_by_sample.get(sample)

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
    out_dir = f"dnaqc_summary/workspace={workspace}/workflow_id={workflow_id}/biosample={sample}"
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   os.path.join(out_dir, "output.parquet"))

    all_summaries.append(summary)

    # Rewrite TSV with final QC_Status, Score, and all columns (MAPD/SKEW from cnvSummarizer)
    score_str = f"{int(composite_score)}" if pd.notna(composite_score) else "NA"
    mapd = cnv_mapd_by_sample.get(sample)
    skew = cnv_skew_by_sample.get(sample)
    mapd_str = f"{mapd:.4f}" if mapd is not None else "NA"
    skew_str = f"{skew:.4f}" if skew is not None else "NA"
    with open(f"{sample}_selected_metrics_mqc.txt", 'w') as fh:
        fh.write("# id: 'dnaqc_summary'\\n")
        fh.write("# plot_type: 'table'\\n")
        fh.write("# section_name: 'QC Summary'\\n")
        fh.write("# description: 'Per-sample alignment, CNV noise, and library complexity QC metrics with pass/fail status.'\\n")
        fh.write("# pconfig:\\n")
        fh.write("#   id: 'dnaqc_summary_table'\\n")
        fh.write("sample_name\\tQC_Status\\tScore\\tpreseq_count\\tPCT_CHIMERAS\\tchrM_pct\\tMAPD_CNV\\tSKEW_CNV\\tTotalReads\\tFinalReads\\tinsert_size\\tfastp_q30_rate\\t%_PF\\t%_Optical_Dups\\t%_PCR_Dups\\n")
        fh.write(f"{sample}\\t{qc_status}\\t{score_str}\\t")
        fh.write(f"{summary.get('preseq_count', 0)}\\t{summary.get('align_pct_chimeras', 0):.4f}\\t")
        fh.write(f"{summary.get('pct_chrm', 0):.4f}\\t{mapd_str}\\t{skew_str}\\t")
        fh.write(f"{summary.get('total_reads', 0)}\\t{summary.get('final_reads', 0)}\\t")
        fh.write(f"{summary.get('insert_size', 0):.0f}\\t{summary.get('fastp_q30_rate', 0):.4f}\\t")
        fh.write(f"{summary.get('pct_pf', 0):.4f}\\t{summary.get('pct_optical_duplicates', 0):.4f}\\t{summary.get('pct_pcr_duplicates', 0):.4f}\\n")

    print(f"  {sample}: QC {qc_status} ({score_str}), CompositeScore={composite_score}")

# Write combined TSV of all samples (for nf-test validation)
pd.DataFrame(all_summaries).to_csv("dnaqc_all_metrics.tsv", sep="\\t", index=False)
pd.DataFrame(all_summaries)[['biosample','qc_status','pipeline','pipeline_version']].rename(
    columns={'biosample':'biosampleName'}).to_csv('per_biosample_status.csv', index=False)
print(f"Written {len(json_files)} parquets with R-derived QC status")
PYEOF
    echo "  ✓ Finished Step 8: parquets written with R consensus QC status"

    # Prepend MultiQC custom-content headers so the SummaryTable renders as a table in the report
    tmp=\$(mktemp)
    printf '# id: "qc_score_distribution"\n# section_name: "Total Usable Cells"\n# description: "Distribution of samples across composite QC score categories. Scores range from 0 (poor quality) to 5 (high quality)."\n# plot_type: "table"\n# pconfig:\n#   id: "qc_score_dist_table"\n#   title: "Total Usable Cells"\n' | cat - DNA-QC_ConsensusScores_SummaryTable_mqc.txt > "\$tmp" && mv "\$tmp" DNA-QC_ConsensusScores_SummaryTable_mqc.txt

    echo ""
    echo "=============================================="
    echo "QC_PLOTS: Successfully completed!"
    echo "=============================================="
    """
}
// Description: Convert SegCopy bin-level CNV data to long-format Parquet
// Output: One parquet per biosample with (chr, start, end, ploidy, raw_count) per bin
// Includes bin_size and read_length for reference selection metadata
// ============================================================================
process GINKGO_BINS_TO_PARQUET {
    tag "cnv_bins_to_parquet"

    input:
    path(segcopy_file)
    path(raw_counts_merged_file)  // Merged binUnsorted.outdata file (matches SegCopy structure)
    val(bin_size)
    val(read_length)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(user)

    output:
    path("cnv_summary/workspace=${workspace}/workflow_id=${workflow_id}/biosample=*/output.parquet"), emit: cnv_bins_parquet

    script:
    """
    python3 - << 'PYEOF'
import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

workspace = "${workspace}"
workflow_id = "${workflow_id}"
dataset_id = "${dataset_id}"
pipeline_version = "${pipeline_version}"
user = "${user}"
bin_size = int("${bin_size}")
read_length = int("${read_length}")

# Read SegCopy (wide format: CHR, START, END, sample1_ploidy, sample2_ploidy, ...)
df_seg = pd.read_csv("${segcopy_file}", sep="\\t")

# Read merged raw counts file (has header with sample names, then raw counts per bin)
df_raw = pd.read_csv("${raw_counts_merged_file}", sep="\\t")

# Get sample columns from SegCopy (everything except CHR, START, END)
sample_cols = [c for c in df_seg.columns if c not in ["CHR", "START", "END"]]

# Get sample columns from raw counts file
raw_sample_cols = list(df_raw.columns)

print(f"SegCopy has {len(df_seg)} bins and {len(sample_cols)} samples: {sample_cols}")
print(f"Raw counts has {len(df_raw)} bins and {len(raw_sample_cols)} samples: {raw_sample_cols}")

# Verify dimensions match
if len(df_seg) != len(df_raw):
    raise ValueError(f"Bin count mismatch: SegCopy has {len(df_seg)} bins, raw counts has {len(df_raw)} bins")
if len(sample_cols) != len(raw_sample_cols):
    raise ValueError(f"Sample count mismatch: SegCopy has {len(sample_cols)} samples, raw counts has {len(raw_sample_cols)} columns")

# Convert to long format and write per-biosample parquet
for i, sample_col in enumerate(sample_cols):
    biosample = sample_col.replace("_sorted", "")
    
    # Find matching raw counts column (sample names should match)
    raw_col = raw_sample_cols[i]
    
    # Create long-format dataframe for this sample (ploidy + raw_count)
    df_long = pd.DataFrame({
        "biosample": biosample,
        "dataset_id": dataset_id,
        "pipeline": "basej-dnaqc",
        "pipeline_version": pipeline_version,
        "user": user,
        "bin_size": bin_size,
        "read_length": read_length,
        "chr": df_seg["CHR"],
        "start": df_seg["START"],
        "end": df_seg["END"],
        "ploidy": df_seg[sample_col],
        "raw_count": df_raw[raw_col]
    })
    
    # Enforce schema types for cnv_summary
    _str_cols    = ['biosample','dataset_id','pipeline','pipeline_version','user','chr']
    _double_cols = ['ploidy']
    _bigint_cols = ['bin_size','read_length','start','end','raw_count']
    for _col in _double_cols:
        if _col in df_long.columns:
            df_long[_col] = pd.to_numeric(df_long[_col], errors='coerce').astype('float64')
    for _col in _bigint_cols:
        if _col in df_long.columns:
            df_long[_col] = pd.to_numeric(df_long[_col], errors='coerce').astype('Int64')
    for _col in _str_cols:
        if _col in df_long.columns:
            df_long[_col] = df_long[_col].astype('string')
    # Write parquet
    out_dir = f"cnv_summary/workspace={workspace}/workflow_id={workflow_id}/biosample={biosample}"
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df_long, preserve_index=False),
                   os.path.join(out_dir, "output.parquet"))

print(f"Created CNV summary parquet files for {len(sample_cols)} samples with ploidy and raw_count")
PYEOF
    """
}

// ============================================================================
// PROCESS: MULTIQC
// Description: Generate QC report combining all metrics
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
    path("multiqc_report_data"), emit: data  // MultiQC v1.28+ uses report name prefix

    script:
    """
    cat > multiqc_config.yaml << EOF
custom_logo_title: 'BioSkryb Genomics'
custom_logo: bioskryb_logo-tagline.png
custom_logo_width: 260

title: "basej-dnaqc v${pipeline_version}"
report_header_info:
  - Dataset ID: "${dataset_id}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false
skip_generalstats: true

# Sample name cleaning
fn_clean_exts:
  - ".fastq.gz"
  - ".fq.gz"
  - ".bam"
  - "_sorted"
  - "_fastp"
  - "_sentieonmetrics"
  - "_trim"
  - ".txt"
  - ".json"

extra_fn_clean_exts:
  - "_R1"
  - "_R2"
  - "_sorted"
  - "_fastp"
  - "_sentieonmetrics"
  - "_trim"
  - ".dedup"

# Module order - custom sections first
module_order:
  - custom_content
  - fastp
  - picard

report_section_order:
  dnaqc_summary:
    order: 1000
  QC_composition_mqc.jpg:
    order: 900
  qc_score_distribution:
    order: 800

table_cond_formatting_rules:
  QC_Status:
    pass:
      - s_eq: "PASS"
    warn:
      - s_eq: "Borderline"
    fail:
      - s_eq: "FAIL"

custom_data:
  QC_composition_mqc.jpg:
    section_name: "QC Composition"
    description: "Per-sample QC metric distributions grouped by cluster. Each panel shows a key metric across all cells."
  dnaqc_summary:
    headers:
      QC_Status:
        title: "QC Status"
        description: "Composite QC status: PASS, Borderline, or FAIL"
        placement: 100
      Score:
        title: "Score"
        description: "Composite QC score (0-5, e.g. 3)"
        placement: 110
      preseq_count:
        title: "Preseq Count"
        description: "Estimated library complexity (number of distinct molecules). Source: preseq."
        format: "{:,.0f}"
        placement: 120
      PCT_CHIMERAS:
        title: "Proportion Chimeras"
        description: "Proportion of chimeric read pairs (0-1). Source: Sentieon AlignmentStat PCT_CHIMERAS."
        format: "{:.4f}"
        placement: 130
      chrM_pct:
        title: "Proportion chrM"
        description: "Proportion of reads mapping to the mitochondrial genome (0-1)."
        format: "{:.4f}"
        placement: 140
      MAPD_CNV:
        title: "MAPD"
        description: "Median Absolute Pairwise Difference of CNV bins — measure of copy-number noise. Lower = better. Source: Ginkgo."
        format: "{:.4f}"
        placement: 150
      SKEW_CNV:
        title: "CNV Skew"
        description: "Skew of CNV bin coverage distribution. Values near 0 indicate balanced coverage. Source: Ginkgo."
        format: "{:.4f}"
        placement: 160
      TotalReads:
        title: "Total Reads"
        description: "Total read pairs before trimming and filtering."
        format: "{:,.0f}"
        placement: 170
      FinalReads:
        title: "Final Reads"
        description: "Read pairs after fastp trimming and deduplication."
        format: "{:,.0f}"
        placement: 180
      insert_size:
        title: "Insert Size (bp)"
        description: "Median insert size in base pairs. Source: Sentieon InsertSizeMetricAlgo."
        format: "{:.0f}"
        placement: 190
      fastp_q30_rate:
        title: "Proportion Q30"
        description: "Proportion of bases with Phred quality >= 30 (0-1). Source: fastp (FASTQ input) or Sentieon QualityYield (CRAM input)."
        format: "{:.4f}"
        placement: 200
EOF

    multiqc . -n multiqc_report.html -c multiqc_config.yaml --force
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
    
    // Auto-detect platform from CSV: parse with splitCsv and check whether
    // *any* row has a non-empty 'cram' column (handles blank first rows and
    // quoted fields that would trip up a naive split(',') approach).
    def csv_file = file(params.input_csv)
    def csv_rows = csv_file.splitCsv(header: true)
    def is_ultima = csv_rows.any { row -> row.cram?.trim() }
    def detected_platform = is_ultima ? 'Ultima' : 'Illumina'

    log.info "Auto-detected platform: ${detected_platform}"
    
    // Parse input CSV based on platform
    if (is_ultima) {
        // Only requires: biosampleName, cram
        // Index file (.cram.crai) is auto-discovered by Nextflow
        ch_reads = channel.fromPath(params.input_csv, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row.cram }
            .map { row ->
                def cram = file(row.cram)
                def crai = row.crai ? file(row.crai) : file(row.cram + '.crai')
                [row.biosampleName, cram, crai]
            }

        ch_reads.view { sample -> "Processing Ultima sample: ${sample[0]}" }
        ch_reads.ifEmpty { exit 1, "ERROR: No CRAM files specified in --input_csv" }

        // CRAM path: subsample to target reads (same as Illumina SEQKIT_SAMPLE)
        SAMTOOLS_SUBSAMPLE_CRAM(
            ch_reads.map { sample_id, cram_file, crai_file -> [sample_id, cram_file, crai_file, params.n_reads] },
            params.samtools_seed,
            file(params.ref_fasta, checkIfExists: true)
        )

        // Pass subsampled CRAMs to alignment
        ch_aligned_reads = SAMTOOLS_SUBSAMPLE_CRAM.out.cram
            .map { sample_id, cram_file, crai_file -> [sample_id, [cram_file, crai_file], "CRAM"] }
        
        // Store read count metrics for QC_PLOTS (matching Illumina metadata)
        ch_fastp_metrics = channel.empty()
        ch_readcount_metrics = SAMTOOLS_SUBSAMPLE_CRAM.out.read_counts_file.collect()
        // Use value channels (never close) so QC_PLOTS process is not terminated early
        ch_fastp_metrics_qc = Channel.value([])
        ch_readcount_metrics_qc = ch_readcount_metrics
        
        // Per-sample read counts for zero-read filtering
        ch_read_counts = SAMTOOLS_SUBSAMPLE_CRAM.out.read_counts
        
    } else {
        // ========== Illumina/Element path: FASTQ input ==========
        // Multi-lane support: read1/read2 can contain pipe-delimited ("|") paths
        //   Single-lane: biosampleName,s3://bucket/R1.fastq.gz,s3://bucket/R2.fastq.gz
        //   Multi-lane:  biosampleName,s3://bucket/L001_R1.fq.gz|s3://bucket/L002_R1.fq.gz,s3://bucket/L001_R2.fq.gz|s3://bucket/L002_R2.fq.gz
        ch_reads_csv = channel.fromPath(params.input_csv, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row.read1 && row.read2 }

        ch_reads_csv.ifEmpty { exit 1, "ERROR: No reads specified in --input_csv" }

        // Branch on pipe character to detect multi-lane vs single-lane
        ch_reads_branched = ch_reads_csv
            .branch {
                multilane: it.read1.contains('|')
                singlelane: true
            }

        // Single-lane: one R1, one R2 — pass through directly
        ch_singlelane = ch_reads_branched.singlelane
            .map { row -> [row.biosampleName, [file(row.read1), file(row.read2)]] }

        ch_singlelane.view { sample -> "Processing ${detected_platform} sample: ${sample[0]}" }

        // Multi-lane: split on "|", collect file objects, flatten for cat
        MERGE_MULTILANE_FASTQ(
            ch_reads_branched.multilane
                .map { row ->
                    def r1_files = row.read1.tokenize('|').collect { file(it.trim()) }
                    def r2_files = row.read2.tokenize('|').collect { file(it.trim()) }
                    [row.biosampleName, r1_files + r2_files]
                }
        )

        ch_fastq_input = ch_singlelane.mix(MERGE_MULTILANE_FASTQ.out.reads)

        // FASTQ path: run SEQKIT and FASTP
        SEQKIT_SAMPLE(
            ch_fastq_input.map { sample_id, reads -> [sample_id, reads, params.n_reads] },
            params.seqkit_sample_seed
        )
        
        FASTP_TRIM(SEQKIT_SAMPLE.out.reads)
        
        ch_aligned_reads = FASTP_TRIM.out.reads
            .map { sample_id, reads -> [sample_id, reads, "FASTQ"] }
        
        // Store metric channels for later use
        ch_fastp_metrics = FASTP_TRIM.out.json_flat.collect()
        ch_readcount_metrics = SEQKIT_SAMPLE.out.read_counts_file.collect()
        ch_fastp_metrics_qc = ch_fastp_metrics
        ch_readcount_metrics_qc = ch_readcount_metrics
        
        // Per-sample read counts for zero-read filtering
        ch_read_counts = SEQKIT_SAMPLE.out.read_counts
    }
    
    // Step 4: Alignment + Dedup + Metrics (handles both FASTQ and CRAM)
    // Choose between Sentieon (proprietary, faster) and open-source path (BWA-MEM2 + samtools markdup + Picard)
    if (params.pipeline_tool == 'sentieon') {
        SENTIEON_ALIGN_DEDUP_METRICS(
            ch_aligned_reads,
            params.reference,
            params.dbsnp,
            params.dbsnp_index,
            params.mills,
            params.mills_index,
            params.onekg_omni,
            params.onekg_omni_index,
            detected_platform,
            params.intervals
        )
        ch_aligned_bam     = SENTIEON_ALIGN_DEDUP_METRICS.out.bam
        ch_align_metrics   = SENTIEON_ALIGN_DEDUP_METRICS.out.metrics
        ch_dedup_metrics   = SENTIEON_ALIGN_DEDUP_METRICS.out.dedup_metrics
        ch_metrics_flat    = SENTIEON_ALIGN_DEDUP_METRICS.out.metrics_flat
    } else {
        // Open-source path: BWA-MEM2 → samtools markdup → Picard metrics.
        // Uses the pre-built BWA-MEM2 index staged at params.bwamem2_reference
        // (built separately by the nf-bwamem2-index-module workflow).
        BWAMEM2_ALIGN(
            ch_aligned_reads,
            params.bwamem2_reference,
            detected_platform
        )
        SAMTOOLS_MARKDUP(BWAMEM2_ALIGN.out.bam)
        PICARD_METRICS(
            SAMTOOLS_MARKDUP.out.bam,
            params.reference,
            params.intervals
        )
        ch_aligned_bam     = SAMTOOLS_MARKDUP.out.bam
        // Combine Picard metrics with samtools dedup metrics into a unified channel
        ch_align_metrics   = PICARD_METRICS.out.metrics
            .join(SAMTOOLS_MARKDUP.out.dedup_metrics)
            .map { sample, picard_files, dedup_file ->
                def all_files = (picard_files instanceof List ? picard_files : [picard_files]) + [dedup_file]
                [sample, all_files]
            }
        ch_dedup_metrics   = SAMTOOLS_MARKDUP.out.dedup_metrics
        ch_metrics_flat    = PICARD_METRICS.out.metrics_flat
            .mix(SAMTOOLS_MARKDUP.out.dedup_metrics_flat)
    }

    // Step 5: Library complexity
    PRESEQ(ch_aligned_bam)
    
    // Step 6: CNV calling with Ginkgo (for supported genomes)
    // Uses GINKO_NOPUBLISH workflow with proper container separation
    ch_ginkgo_rds = channel.empty()
    ch_ginkgo_segcopy = channel.empty()
    ch_ginkgo_metrics = channel.empty()
    ch_ginkgo_raw_counts = channel.empty()
    ch_ginkgo_cnv_plots = channel.empty()
    ch_cnv_summary_parquet = channel.empty()
    
    if (params.genome in ["GRCh38", "GRCm39", "ARSUCD2"]) {
        // Filter out zero-read samples before Ginkgo CNV
        // Zero-read BAMs produce degenerate RDS/SegCopy that crash R scripts
        // Samples with 0 reads will still appear in QC_PLOTS/parquet/MultiQC with NULL CNV values
        ch_bam_for_ginkgo = ch_aligned_bam
            .join(ch_read_counts)   // [sample, bam, bai, total_reads, final_reads]
            .filter { sample, bam, bai, total_reads, final_reads ->
                if (total_reads.toString().toLong() == 0) {
                    log.warn "Skipping Ginkgo CNV for sample '${sample}' — 0 reads"
                    return false
                }
                return true
            }
            .map { sample, bam, bai, total_reads, final_reads -> [sample, bam, bai] }

        // Setup Ginkgo references
        // Use channel.of(file()) instead of channel.fromPath() to avoid S3 listing at startup
        // channel.fromPath() actively lists S3 before credentials are initialized, returning empty channels
        // channel.of(file()) creates lazy path objects; S3 access is deferred to task staging time
        ch_binref = channel.of(file(params.ginko_ref_dir + "variable_" + params.bin_size + "_" + params.read_length + "_bwa"))
        ch_gcref = channel.of(file(params.ginko_ref_dir + "GC_variable_" + params.bin_size + "_" + params.read_length + "_bwa"))
        ch_boundsref = channel.of(file(params.ginko_ref_dir + "bounds_variable_" + params.bin_size + "_" + params.read_length + "_bwa"))
        
        // Run Ginkgo CNV workflow (containers defined inline in nextflow.config)
        GINKO_NOPUBLISH(
            ch_bam_for_ginkgo,
            params.bin_size,
            ch_binref,
            ch_gcref,
            ch_boundsref,
            params.min_ploidy,
            params.max_ploidy,
            params.min_bin_width,
            params.is_haplotype
        )
        
        ch_ginkgo_rds = GINKO_NOPUBLISH.out.rds
        ch_ginkgo_segcopy = GINKO_NOPUBLISH.out.segcopy
        ch_ginkgo_metrics = GINKO_NOPUBLISH.out.metrics
        ch_ginkgo_raw_counts = GINKO_NOPUBLISH.out.raw_counts
        ch_ginkgo_cnv_plots = GINKO_NOPUBLISH.out.graph
        
        // Convert bin-level CNV data to long-format Parquet (per biosample)
        GINKGO_BINS_TO_PARQUET(
            GINKO_NOPUBLISH.out.segcopy,
            GINKO_NOPUBLISH.out.raw_counts_merged,
            params.bin_size,
            params.read_length,
            params.dataset_id,
            params.workspace,
            params.workflow_id,
            workflow.manifest.version,
            params.pipeline_user
        )
        ch_cnv_summary_parquet = GINKGO_BINS_TO_PARQUET.out.cnv_bins_parquet
    }
    
    // CNV is required - fail if outputs are missing
    ch_ginkgo_metrics_file = ch_ginkgo_metrics
        .ifEmpty { exit 1, "ERROR: Ginkgo CNV workflow failed - CNV metrics are required for DNA-QC pipeline" }
    ch_ginkgo_segcopy_file = ch_ginkgo_segcopy
        .ifEmpty { exit 1, "ERROR: Ginkgo CNV workflow failed - CNV segcopy is required for DNA-QC pipeline" }
    
    // Step 7: Consolidated QC_PLOTS - parses metrics, merges, and generates plots
    // ch_fastp_metrics and ch_readcount_metrics are created based on platform
    // QC_PLOTS script uses glob patterns to find files, works with empty inputs
    QC_PLOTS(
        ch_ginkgo_rds,
        ch_ginkgo_segcopy_file,
        ch_ginkgo_metrics_file,
        file(params.input_csv),
        file("${projectDir}/assets/plot_qc_config.json"),
        ch_metrics_flat.collect(),
        ch_fastp_metrics_qc,
        PRESEQ.out.complexity_flat.collect(),
        ch_readcount_metrics_qc,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        params.pipeline_user
    )
    
    // Step 8: MultiQC Report
    // Mix all metric sources as flat channels, then collect into single list for MULTIQC
    // Empty channels (from Ultima) naturally produce no items, which is fine
    ch_mqc_inputs = ch_metrics_flat
        .mix(ch_fastp_metrics)
        .mix(PRESEQ.out.complexity_flat)
        .mix(QC_PLOTS.out.mqc_metrics)
        .mix(QC_PLOTS.out.composition_jpg)
        .mix(QC_PLOTS.out.cnv_quadrants_jpg)
        .mix(QC_PLOTS.out.consensus_summary)
        .collect()
    
    MULTIQC(
        ch_mqc_inputs,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        file("${projectDir}/assets/bioskryb_logo-tagline.png")
    )

    publish:
    bam_files = ch_aligned_bam
        .map { sample_name, bam, bai -> [biosampleName: sample_name, bam: bam, bai: bai] }
    
    metrics_files = ch_align_metrics
        .map { sample_name, metrics -> [biosampleName: sample_name, metrics: metrics] }
    
    dedup_metrics = ch_dedup_metrics
        .map { sample_name, metrics -> [biosampleName: sample_name, metrics: metrics] }
    
    preseq_complexity = PRESEQ.out.complexity
        .map { sample_name, preseq -> [biosampleName: sample_name, preseq: preseq] }
    
    // Single unified dnaqc_summary with ALL metrics (from consolidated QC_PLOTS)
    // Contains: Sentieon + Fastp + Preseq + Ginkgo CNV summary per biosample
    dnaqc_summary = QC_PLOTS.out.parquet
    
    // Combined QC metrics for all samples (from consolidated QC_PLOTS)
    qc_metrics_combined = QC_PLOTS.out.combined_metrics

    // Full-schema TSV of all samples for nf-test validation
    dnaqc_summary_tsv = QC_PLOTS.out.summary_tsv
    
    // QC_PLOTS outputs - all consolidated
    qc_plots = QC_PLOTS.out.allmetrics_with_cnv
        .mix(QC_PLOTS.out.consensus_scores)
        .mix(QC_PLOTS.out.consensus_summary)
        .mix(QC_PLOTS.out.consensus_group_summary)
        .mix(QC_PLOTS.out.composition_pdf)
        .mix(QC_PLOTS.out.cnv_quadrants_pdf)
        .mix(QC_PLOTS.out.composition_jpg)
        .mix(QC_PLOTS.out.cnv_quadrants_jpg)
        .collect()
    
    // CNV-Summary: Ploidy calls per biosample per bin (long format)
    cnv_summary = ch_cnv_summary_parquet
    
    ginkgo_rds = ch_ginkgo_rds
    ginkgo_segcopy = ch_ginkgo_segcopy
    cnv_plots_per_cell = ch_ginkgo_cnv_plots
    per_biosample_status = QC_PLOTS.out.per_biosample_status
    multiqc_report = MULTIQC.out.report
}

// ============================================================================
// OUTPUT CONFIGURATION
// ============================================================================
output {
    bam_files {
        path "bam/${params.workspace}/dna/tool=${params.pipeline_tool == 'sentieon' ? 'sentieon-202503-02' : 'bwa-mem2'}/pipeline=dnaqc"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/bam.csv"
            header true
        }
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "bam",
             tool: "${params.pipeline_tool == 'sentieon' ? 'sentieon' : 'bwa-mem2'}".toString(),
             reference: params.genome
    }
    
    metrics_files {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/sentieon"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/metrics.csv"
            header true
        }
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "qc_metrics",
             tool: "${params.pipeline_tool == 'sentieon' ? 'sentieon' : 'picard'}".toString(),
             reference: params.genome
    }
    
    dedup_metrics {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/dedup"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/dedup_metrics.csv"
            header true
        }
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "dedup_metrics",
             tool: "${params.pipeline_tool == 'sentieon' ? 'sentieon' : 'samtools'}".toString(),
             reference: params.genome
    }
    
    preseq_complexity {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/preseq"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "preseq",
             tool: "preseq"
    }

    // Unified dnaqc_summary parquet - contains ALL metrics per biosample
    // (Sentieon, Fastp, Preseq, and Ginkgo CNV summary when available)
    dnaqc_summary {
        path "tables"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "dnaqc_summary",
             reference: params.genome
    }

    // Combined QC metrics TSV for all samples (used by QC_PLOTS)
    qc_metrics_combined {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/qc_metrics"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "qc_metrics_combined"
    }

    // Full-schema TSV mirroring the parquet schema (for nf-test validation)
    dnaqc_summary_tsv {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/qc_metrics"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "dnaqc_summary_tsv"
    }

    // QC_PLOTS outputs - all consolidated
    qc_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/qc_plots"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "qc_plots"
    }

    // CNV-Summary: Ploidy calls per biosample per bin (long format)
    // Includes bin_size and read_length for reference selection metadata
    cnv_summary {
        path "tables"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_summary",
             tool: "ginkgo",
             reference: params.genome
    }

    ginkgo_rds {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/cnv/ginkgo_rds"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_rds",
             tool: "ginkgo"
    }

    ginkgo_segcopy {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/cnv/ginkgo_segcopy"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_segcopy",
             tool: "ginkgo"
    }
    
    // Per-cell CNV profile plots (tar.gz archive of one JPEG per biosample)
    cnv_plots_per_cell {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/cnv/ginkgo_cnv_plots"
        tags workspace: params.workspace,
             dataset_id: params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline: workflow.manifest.name,
             molecule_type: "dna",
             artifact: "cnv_plots",
             tool: "ginkgo"
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
