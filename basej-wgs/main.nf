nextflow.enable.dsl=2

// ============================================================================
// BASEJ-WGS: Single-Cell WGS/WES QC Pipeline
// ============================================================================
// Description: Alignment, deduplication, BQSR, and QC metrics collection.
//              Supports both WGS (whole genome) and WES (whole exome) modes
//              via params.mode = 'wgs' | 'exome'.
// Outputs: Per-biosample Parquet files for Athena queries + MultiQC report
// ============================================================================

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
// Description: Optional subsampling — count reads with zcat + wc,
//              then if total reads > max_total_reads, seqkit sample -p (no two-pass -2) caps reads.
//              Paired: R1 then R2 sequentially with same PROPORTION/seed; -j 2.
//              skip_subsampling=true by default.
// ============================================================================
process SEQKIT_SAMPLE {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(max_total_reads)
    val(seqkit_sample_seed)

    output:
    tuple val(sample_name), path("${sample_name}_subsampled_R*.fastq.gz"), emit: reads

    script:
    def r1 = reads[0]
    def r2 = reads.size() == 2 ? reads[1] : ""
    def paired = reads.size() == 2
    def count_r2 = paired ? "( zcat '${r2}' | wc -l | awk '{printf \"%.0f\", \$1/4}' > read2.txt ) &" : ""
    def copy_r2 = paired ? "cp '${r2}' '${sample_name}_subsampled_R2.fastq.gz' &" : ""
    def sample_r2 = paired ? "seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j 2 -o '${sample_name}_subsampled_R2.fastq.gz' '${r2}'" : ""
    """
    set -euo pipefail

    # Count reads per mate in parallel (zcat + wc; reads = lines/4)
    ( zcat '${r1}' | wc -l | awk '{printf "%.0f", \$1/4}' > read1.txt ) &
    ${count_r2}
    wait

    R1=\$(cat read1.txt)
    if [ "${paired}" = "true" ]; then
      R2=\$(cat read2.txt)
      TOTAL_READS=\$((R1 + R2))
    else
      TOTAL_READS=\$R1
    fi
    echo "Total reads (all mates): \$TOTAL_READS"

    TARGET=${max_total_reads}

    if [ "\$TOTAL_READS" -le "\$TARGET" ]; then
        echo "No subsampling (total <= \$TARGET); copying inputs..."
        cp '${r1}' '${sample_name}_subsampled_R1.fastq.gz' &
        ${copy_r2}
        wait
    else
        echo "Subsampling to \$TARGET total reads with seqkit sample -p -j 2 (R1 then R2, no two-pass)..."
        export PROPORTION=\$(awk -v t="\$TARGET" -v tot="\$TOTAL_READS" 'BEGIN { printf "%.18f", t/tot }')
        seqkit sample -p \$PROPORTION -s ${seqkit_sample_seed} -j 2 \\
          -o '${sample_name}_subsampled_R1.fastq.gz' '${r1}'
        ${paired ? sample_r2 : ""}
        if [ "${paired}" = "true" ]; then
          test -f '${sample_name}_subsampled_R1.fastq.gz' && test -f '${sample_name}_subsampled_R2.fastq.gz' || { echo "ERROR: subsampling did not produce both mates"; exit 1; }
        fi
    fi
    """
}

// ============================================================================
// PROCESS: SAMTOOLS_SUBSAMPLE_CRAM
// Description: Optional subsampling of CRAM files for Ultima platform —
//              count reads with samtools view -c, then if total > max_total_reads,
//              subsample with samtools view -s (deterministic seeding).
//              Output converges to reads channel with read_type="CRAM" tag.
// ============================================================================
process SAMTOOLS_SUBSAMPLE_CRAM {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(cram), path(crai), val(max_total_reads)
    val(samtools_seed)

    output:
    tuple val(sample_name), path("${sample_name}.cram"), path("${sample_name}.cram.crai"), emit: reads
    path("${sample_name}_read_counts.txt"), emit: read_counts_file

    script:
    """
    # Count total reads in CRAM
    export TOTAL_READS=\$(samtools view -c '${cram}')
    echo "Total reads: \$TOTAL_READS"
    
    # Calculate subsample proportion (same logic as SEQTK_SAMPLE)
    PROPORTION=\$(awk -v total=\$TOTAL_READS -v target=${max_total_reads} 'BEGIN { p=target/total; if(p>1) p=1; print p }')
    echo "Subsample proportion: \$PROPORTION"
    
    if [ \$(awk -v p="\$PROPORTION" 'BEGIN{if (p >= 1) print 1; else print 0}') -eq 1 ]; then
        echo "No subsampling needed, copying CRAM..."
        if [ '${cram}' != '${sample_name}.cram' ]; then
            cp '${cram}' '${sample_name}.cram'
            cp '${crai}' '${sample_name}.cram.crai'
        fi
        export FINAL_READS=\$TOTAL_READS
    else
        echo "Subsampling to ${max_total_reads} reads..."
        # Use samtools view -s for deterministic subsampling with seed.
        # Pin CRAM 3.0: samtools >=1.22 defaults to CRAM 3.1, whose codecs
        # the DeepVariant 1.8.0 htslib cannot decode (Failure to decode slice).
        samtools view -s ${samtools_seed}.\$PROPORTION -C --output-fmt cram,version=3.0 -o '${sample_name}.cram' '${cram}'
        samtools index '${sample_name}.cram' '${sample_name}.cram.crai'
        # Verify final read count
        export FINAL_READS=\$(samtools view -c '${sample_name}.cram')
    fi
    
    echo "Final reads: \$FINAL_READS"
    
    # Write read counts to file for QC_PLOTS (matching SEQTK_SAMPLE format)
    echo "\$TOTAL_READS" > '${sample_name}_read_counts.txt'
    echo "\$FINAL_READS" >> '${sample_name}_read_counts.txt'
    """
}

// ============================================================================
// PROCESS: SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION
// Description: Coverage-targeted CRAM subsampling for Ultima fan-out runs.
//              Unlike SAMTOOLS_SUBSAMPLE_CRAM (which counts reads first to
//              hit a max_total_reads target), this process is given a
//              pre-computed proportion derived in the workflow from
//              `target_coverage / mean_coverage` (mean_coverage is supplied
//              per-row in the input CSV). When a biosample's native coverage
//              is at or below the target, proportion is 1.0 and the CRAM is
//              passed through unchanged (no upsampling) and labeled with its
//              actual floored coverage. Output CRAM is named
//              `{fanout_name}.cram` where fanout_name encodes the coverage
//              (e.g. `HG002_pta_001_cov10x`). Dedup flags (0x400)
//              are preserved by `samtools view -s` so downstream metrics
//              still see the upstream Ultima `demux` duplicate marks.
// ============================================================================
process SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION {
    tag "${fanout_name}"

    input:
    tuple val(fanout_name), path(cram), path(crai), val(proportion)
    val(samtools_seed)

    output:
    tuple val(fanout_name), path("${fanout_name}.cram"), path("${fanout_name}.cram.crai"), emit: reads

    script:
    """
    set -euo pipefail
    # Pin CRAM 3.0: samtools >=1.22 defaults to CRAM 3.1, whose codecs the
    # DeepVariant 1.8.0 htslib cannot decode (Failure to decode slice).
    # Sentieon's newer htslib reads 3.1 fine, which is why wgsqc metrics
    # succeed on the same files while DeepVariant make_examples fails.
    if awk -v p="${proportion}" 'BEGIN { exit !(p >= 1.0) }'; then
        # proportion >= 1.0: native coverage is at/below target, cannot
        # upsample. Pass through unchanged (only transcode to CRAM 3.0).
        echo "Pass-through ${cram} at native coverage (proportion ${proportion} >= 1.0, no subsampling)"
        samtools view -C --output-fmt cram,version=3.0 -@ ${task.cpus} \\
            -o '${fanout_name}.cram' '${cram}'
    else
        # samtools view -s expects SEED.FRAC where FRAC is the fractional digits
        # of the proportion. Format proportion to 6 decimal places, strip leading
        # "0.".
        FRAC=\$(awk -v p="${proportion}" 'BEGIN { printf "%.6f", p }' | cut -d. -f2)
        echo "Subsampling ${cram} to proportion ${proportion} (samtools -s ${samtools_seed}.\$FRAC)"
        samtools view -s ${samtools_seed}.\$FRAC -C --output-fmt cram,version=3.0 -@ ${task.cpus} \\
            -o '${fanout_name}.cram' '${cram}'
    fi
    samtools index -@ ${task.cpus} '${fanout_name}.cram' '${fanout_name}.cram.crai'

    FINAL_READS=\$(samtools view -c -@ ${task.cpus} '${fanout_name}.cram')
    echo "Final ${fanout_name}: \$FINAL_READS reads"
    """
}

// ============================================================================
// PROCESS: SENTIEON_ALIGN_DEDUP
// Description: BWA-MEM alignment + deduplication (FASTQ path only).
//              Pre-aligned CRAMs (e.g. Ultima) skip this entirely and run
//              SENTIEON_METRICS_CRAM directly — Ultima `demux` already marks
//              duplicates upstream, so re-running Sentieon Dedup is redundant
//              and adds 2+ hours per sample.
// ============================================================================
process SENTIEON_ALIGN_DEDUP {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(read_type)
    path fasta_ref
    val(platform)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("${sample_name}.dedup_sentieonmetrics.txt"), emit: dedup_metrics
    path("${sample_name}.dedup_sentieonmetrics.txt"), emit: dedup_metrics_flat

    script:
    // FASTQ path: reads contains R1 and/or R2 files, proceed with alignment + dedup
    def r1 = reads[0]
    def r2 = reads.size() == 2 ? reads[1] : ""
    """
    set +u
    export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

    export bwt_max_mem=\$([ ${task.memory.toGiga()} -gt 30 ] && echo "30G" || echo "${task.memory.toGiga()}G")

    # BWA-MEM alignment + sort
    sentieon bwa mem -M -Y -K 2500000000 \\
        -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:${platform}" \\
        -t ${task.cpus} '${fasta_ref}/genome.fa' '${r1}' '${r2}' | \\
        sentieon util sort -r '${fasta_ref}/genome.fa' -o '${sample_name}_sorted.bam' -t ${task.cpus} --sam2bam -i -

    samtools index -@ ${task.cpus} '${sample_name}_sorted.bam'

    # LocusCollector for dedup
    sentieon driver -t ${task.cpus} -r '${fasta_ref}/genome.fa' -i ${sample_name}_sorted.bam \\
        --algo LocusCollector --fun score_info '${sample_name}.locuscollector_score.gz'

    sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${sample_name}_sorted.bam \\
        --algo Dedup --score_info ${sample_name}.locuscollector_score.gz \\
        --metrics ${sample_name}.dedup_sentieonmetrics.txt ${sample_name}.bam

    samtools index -@ ${task.cpus} '${sample_name}.bam'

    rm -f '${sample_name}_sorted.bam' '${sample_name}_sorted.bam.bai' \\
          '${sample_name}.locuscollector_score.gz'
    """
}

// ============================================================================
// PROCESS: SENTIEON_DRIVER_METRICS
// Description: Collect alignment, GC bias, insert size, and coverage metrics.
//              WGS mode: WgsMetricsAlgo + CoverageMetrics on full genome intervals
//              Exome mode: CoverageMetrics on target intervals only
// ============================================================================
process SENTIEON_DRIVER_METRICS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path fasta_ref
    path base_metrics_intervals
    path wgs_or_target_intervals
    val mode

    output:
    tuple val(sample_name), path("*sentieonmetrics*"), emit: metrics_tuple
    path("*sentieonmetrics*"), emit: metrics_flat

    script:
    if (mode == 'exome') {
        """
        set +u
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${bam} \\
            --interval ${wgs_or_target_intervals} \\
            --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \\
            --algo AlignmentStat ${sample_name}.alignmentstat_sentieonmetrics.txt \\
            --algo InsertSizeMetricAlgo ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \\
            --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \\
            --algo QualityYield ${sample_name}.qualityyield_sentieonmetrics.txt \\
            --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics

        # WGS metrics not applicable in exome mode
        touch ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt

        export SENTIEON_VER="202308.01"
        echo Sentieon: \$SENTIEON_VER > sentieon_driver_metrics_version.yml
        """
    } else {
        """
        set +u
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${bam} \\
            --interval ${base_metrics_intervals} \\
            --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \\
            --algo AlignmentStat ${sample_name}.alignmentstat_sentieonmetrics.txt \\
            --algo InsertSizeMetricAlgo ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt \\
            --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \\
            --algo QualityYield ${sample_name}.qualityyield_sentieonmetrics.txt \\
            --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics --omit_base_output --omit_locus_stat --omit_sample_stat

        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i ${bam} \\
            --interval ${wgs_or_target_intervals} \\
            --algo WgsMetricsAlgo ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt

        # HS metrics not applicable in WGS mode
        touch ${sample_name}.hsmetricalgo.sentieonmetrics.txt

        export SENTIEON_VER="202308.01"
        echo Sentieon: \$SENTIEON_VER > sentieon_driver_metrics_version.yml
        """
    }
}

// ============================================================================
// PROCESS: SENTIEON_METRICS_CRAM
// Description: Metrics-only path for pre-aligned CRAM input (e.g. Ultima).
//              Runs Sentieon driver directly on the CRAM with no dedup re-run.
//              Ultima `demux` pre-marks duplicates upstream, so WgsMetricsAlgo
//              and CoverageMetrics still honor the existing 0x400 flag.
//              A header check for known dedup @PG entries logs a warning if
//              upstream marking is missing.
//              WGS mode: AlignmentStat + GCBias + MeanQualityByCycle +
//                        WgsMetricsAlgo --include_unpaired + QualityYield
//              Exome mode: AlignmentStat + GCBias + MeanQualityByCycle +
//                          CoverageMetrics + QualityYield
//              Empty placeholder files match the FASTQ-path glob pattern so
//              the existing parser does not need a separate code path.
// ============================================================================
process SENTIEON_METRICS_CRAM {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(cram), path(crai)
    path fasta_ref
    path base_metrics_intervals
    path wgs_or_target_intervals
    val mode

    output:
    tuple val(sample_name), path("*sentieonmetrics*"), emit: metrics_tuple
    path("*sentieonmetrics*"), emit: metrics_flat

    script:
    if (mode == 'exome') {
        """
        set +u
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

        # Compute pct_duplication from the existing 0x400 flag (no re-dedup).
        # Picard's PERCENT_DUPLICATION for unpaired data ≈ primary duplicates / primary mapped.
        samtools flagstat -@ ${task.cpus} '${cram}' > ${sample_name}.flagstat.sentieonmetrics.txt

        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i '${cram}' \\
            --interval ${wgs_or_target_intervals} \\
            --algo AlignmentStat ${sample_name}.alignmentstat_sentieonmetrics.txt \\
            --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \\
            --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \\
            --algo QualityYield ${sample_name}.qualityyield_sentieonmetrics.txt \\
            --algo CoverageMetrics ${sample_name}.cov_sentieonmetrics

        # Placeholders so the parser's glob picks up consistent file lists.
        # InsertSize is paired-only; WgsMetrics is WGS-only.
        # hsmetricalgo is provided by PICARD_COLLECTHSMETRICS in exome mode.
        touch ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt
        touch ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt
        touch ${sample_name}.dedup_sentieonmetrics.txt

        export SENTIEON_VER="202308.01"
        echo Sentieon: \$SENTIEON_VER > sentieon_driver_metrics_version.yml
        """
    } else {
        """
        set +u
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER

        # Compute pct_duplication from the existing 0x400 flag (no re-dedup).
        # Picard's PERCENT_DUPLICATION for unpaired data ≈ primary duplicates / primary mapped.
        samtools flagstat -@ ${task.cpus} '${cram}' > ${sample_name}.flagstat.sentieonmetrics.txt

        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i '${cram}' \\
            --interval ${base_metrics_intervals} \\
            --algo AlignmentStat ${sample_name}.alignmentstat_sentieonmetrics.txt \\
            --algo GCBias --summary ${sample_name}.gcbias_summary.sentieonmetrics.txt ${sample_name}.gcbias.sentieonmetrics.txt \\
            --algo MeanQualityByCycle ${sample_name}.meanqualitybycycle.sentieonmetrics.txt \\
            --algo QualityYield ${sample_name}.qualityyield_sentieonmetrics.txt

        sentieon driver -t ${task.cpus} -r ${fasta_ref}/genome.fa -i '${cram}' \\
            --interval ${wgs_or_target_intervals} \\
            --algo WgsMetricsAlgo --include_unpaired true ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt

        # Placeholders so the parser's glob picks up consistent file lists.
        touch ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt
        touch ${sample_name}.dedup_sentieonmetrics.txt
        touch ${sample_name}.hsmetricalgo.sentieonmetrics.txt
        touch ${sample_name}.cov_sentieonmetrics

        export SENTIEON_VER="202308.01"
        echo Sentieon: \$SENTIEON_VER > sentieon_driver_metrics_version.yml
        """
    }
}

// ============================================================================
// PROCESS: PICARD_COLLECTHSMETRICS
// Description: Collect hybrid selection (exome) metrics - exome mode only
// ============================================================================
process PICARD_COLLECTHSMETRICS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path fasta_ref
    path intervals

    output:
    tuple val(sample_name), path("${sample_name}.hsmetricalgo.sentieonmetrics.txt"), emit: metrics_tuple
    path("${sample_name}.hsmetricalgo.sentieonmetrics.txt"), emit: metrics_flat

    script:
    def avail_mem = task.memory ? task.memory.giga : 3
    """
    picard \\
        -Xmx${avail_mem}g \\
        CollectHsMetrics \\
        --INPUT ${bam} \\
        --OUTPUT ${sample_name}.hsmetricalgo.sentieonmetrics.txt \\
        --BAIT_INTERVALS ${intervals} \\
        --TARGET_INTERVALS ${intervals} \\
        --REFERENCE_SEQUENCE ${fasta_ref}/genome.fa
    """
}

// ============================================================================
// PROCESS: BWAMEM2_ALIGN_DEDUP_METRICS  (open-source, single consolidated process)
// Description: Open-source replacement for SENTIEON_ALIGN_DEDUP + SENTIEON_DRIVER_METRICS.
//              Runs alignment, deduplication AND metrics collection in ONE process so the
//              large WGS BAM never has to be staged/transferred between tasks.
//                1. bwa-mem2 mem  → samtools name-sort            (≈ sentieon bwa + util sort)
//                2. samtools fixmate → sort → markdup (+ Picard DuplicationMetrics) (≈ LocusCollector + Dedup)
//                3. Picard CollectMultipleMetrics (single pass) + CollectWgsMetrics (≈ sentieon driver --algo ...)
//              Output filenames match the Sentieon *sentieonmetrics* glob so the downstream
//              parser (WGS_QC_METRICS_TO_PARQUET) is unchanged.
//              NOTE: bwa-mem2 has no production ARM build, so this process is x86-only.
// ============================================================================
process BWAMEM2_ALIGN_DEDUP_METRICS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads), val(read_type)
    path bwamem2_index_dir
    path fasta_ref
    path base_metrics_intervals
    path wgs_or_target_intervals
    val mode
    val platform

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("${sample_name}.dedup_sentieonmetrics.txt"), emit: dedup_metrics
    path("${sample_name}.dedup_sentieonmetrics.txt"), emit: dedup_metrics_flat
    tuple val(sample_name), path("*sentieonmetrics*"), emit: metrics_tuple
    path("*sentieonmetrics*"), emit: metrics_flat

    script:
    def r1 = reads[0]
    def r2 = reads.size() == 2 ? reads[1] : ""
    def avail_mem = task.memory ? task.memory.giga : 8
    """
    set +u

    # ============ 1. bwa-mem2 alignment → name-sort (for fixmate) ============
    # Pipe alignments straight into a name sort. The previous coordinate sort +
    # index of an intermediate BAM was redundant work: markdup needs a fixmate
    # pass first, and fixmate requires name-grouped input. Uncompressed (-u)
    # intermediates skip pointless compress/decompress between samtools steps.
    BWA_INDEX="${bwamem2_index_dir}/genome.fa"
    bwa-mem2.avx2 mem -M -Y -K 2500000000 \\
        -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:${platform}" \\
        -t ${task.cpus} \$BWA_INDEX '${r1}' '${r2}' | \\
        samtools sort -n -u -@ ${task.cpus} -o ${sample_name}_namesorted.bam -

    # ============ 2. samtools markdup → Picard-style DuplicationMetrics ============
    # fixmate (name-grouped) → coordinate sort → markdup. Same tool/algorithm as
    # before; only the redundant leading coordinate sort + index were dropped.
    samtools fixmate -m -u -@ ${task.cpus} ${sample_name}_namesorted.bam ${sample_name}_fixmate.bam
    samtools sort -u -@ ${task.cpus} -o ${sample_name}_resorted.bam ${sample_name}_fixmate.bam
    # Remove duplicates (-r) to match Sentieon's `Dedup --rmdup` behavior, so the
    # output BAM (and therefore coverage and downstream metrics) is consistent
    # between the open-source and Sentieon paths. `-s` still reports the duplicate
    # counts (computed before removal) used to build the dedup metrics below.
    samtools markdup -r -s -f ${sample_name}.markdup_stats.txt -@ ${task.cpus} \\
        ${sample_name}_resorted.bam ${sample_name}.bam
    samtools index -@ ${task.cpus} ${sample_name}.bam

    # Convert samtools markdup stats → Picard DuplicationMetrics so the parser
    # (READ_PAIRS_EXAMINED / READ_PAIR_DUPLICATES / PERCENT_DUPLICATION) is unchanged.
    # markdup counts reads; Picard counts pairs, so divide paired counts by 2.
    get_stat() { grep -E "^\$1:" ${sample_name}.markdup_stats.txt | head -1 | awk -F': ' '{print \$2+0}'; }
    PAIRED_READS=\$(get_stat "PAIRED")
    DUP_PAIR_READS=\$(get_stat "DUPLICATE PAIR")
    DUP_PAIR_OPT_READS=\$(get_stat "DUPLICATE PAIR OPTICAL")
    UNPAIRED_READS=\$(get_stat "SINGLE")
    UNPAIRED_DUPS=\$(get_stat "DUPLICATE SINGLE")
    UNMAPPED=\$(get_stat "EXCLUDED")
    READ_PAIRS_EXAMINED=\$(awk -v r=\$PAIRED_READS 'BEGIN{printf "%d", r/2}')
    READ_PAIR_DUPLICATES=\$(awk -v r=\$DUP_PAIR_READS 'BEGIN{printf "%d", r/2}')
    READ_PAIR_OPTICAL_DUPLICATES=\$(awk -v r=\$DUP_PAIR_OPT_READS 'BEGIN{printf "%d", r/2}')
    PCT_DUP=\$(awk -v dp=\$READ_PAIR_DUPLICATES -v du=\$UNPAIRED_DUPS -v ep=\$READ_PAIRS_EXAMINED -v eu=\$UNPAIRED_READS \\
        'BEGIN { denom = ep + eu; if (denom > 0) printf "%.6f", (dp + du) / denom; else printf "0" }')

    cat > ${sample_name}.dedup_sentieonmetrics.txt << EOF
## METRICS CLASS	picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
${sample_name}	\${UNPAIRED_READS}	\${READ_PAIRS_EXAMINED}	0	\${UNMAPPED}	\${UNPAIRED_DUPS}	\${READ_PAIR_DUPLICATES}	\${READ_PAIR_OPTICAL_DUPLICATES}	\${PCT_DUP}	0
EOF

    # ============ 3. Picard metric suite — single pass via CollectMultipleMetrics ============
    # Runs AlignmentSummary, InsertSize, GcBias, MeanQualityByCycle and
    # QualityYield as PROGRAM modules in ONE pass over the dedup BAM instead of
    # five separate full-BAM passes. Same Picard collectors → identical metric
    # values; only file names differ, so we rename to the *sentieonmetrics*
    # names the parser and MultiQC expect (MultiQC detects by the internal
    # "## METRICS CLASS picard..." header, not the filename).
    picard -Xmx${avail_mem}g CollectMultipleMetrics \\
        R=${fasta_ref}/genome.fa \\
        I=${sample_name}.bam \\
        O=${sample_name}.multimetrics \\
        PROGRAM=null \\
        PROGRAM=CollectAlignmentSummaryMetrics \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=MeanQualityByCycle \\
        PROGRAM=CollectQualityYieldMetrics

    # Rename single-pass outputs to the expected *sentieonmetrics* file names.
    # Some PROGRAM modules emit no file when there is no eligible data (e.g.
    # CollectInsertSizeMetrics on a sample with zero mapped pairs), so fall back
    # to an empty placeholder — the parser already skips empty metric files.
    mv_or_touch() { if [ -f "\$1" ]; then mv "\$1" "\$2"; else touch "\$2"; fi; }
    mv_or_touch ${sample_name}.multimetrics.insert_size_metrics      ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.gc_bias.summary_metrics  ${sample_name}.gcbias_summary.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.gc_bias.detail_metrics   ${sample_name}.gcbias.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.quality_by_cycle_metrics ${sample_name}.meanqualitybycycle.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.quality_yield_metrics    ${sample_name}.qualityyield_sentieonmetrics.txt

    # AlignmentSummaryMetrics: drop Picard's 3 trailing columns (SAMPLE,
    # LIBRARY, READ_GROUP) so the CATEGORY data rows and header have matching
    # field counts in the parser. Tolerates a missing file (empty input).
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
                of.write("\\t".join(parts[:-3]) + "\\n")
            else:
                of.write(line)
PYEOF
    rm -f ${sample_name}.multimetrics.alignment_summary_metrics

    # ============ 4. WGS coverage metrics (Picard CollectWgsMetrics) ============
    # CollectWgsMetrics is the sole source of every coverage field the parser
    # reads (MEAN_COVERAGE / PCT_*X / PCT_EXC_*). The previous
    # `samtools depth -a | awk` interval summary was never parsed downstream
    # (nor meaningfully consumed by MultiQC), so it is removed — eliminating a
    # full-genome single-threaded depth pass. Quality/coverage cutoffs are
    # pinned to match Sentieon WgsMetricsAlgo defaults so concordance can't drift.
    if [ "${mode}" = "exome" ]; then
        # WgsMetrics N/A in exome mode (HsMetrics provides target coverage).
        touch ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt
    else
        picard -Xmx${avail_mem}g CollectWgsMetrics \\
            R=${fasta_ref}/genome.fa \\
            I=${sample_name}.bam \\
            O=${sample_name}.wgsmetricsalgo.sentieonmetrics.txt \\
            INTERVALS=${wgs_or_target_intervals} \\
            MINIMUM_MAPPING_QUALITY=20 \\
            MINIMUM_BASE_QUALITY=20 \\
            COVERAGE_CAP=250
        # hsmetricalgo N/A in WGS mode
        touch ${sample_name}.hsmetricalgo.sentieonmetrics.txt
    fi

    # Cleanup intermediates
    rm -f ${sample_name}_namesorted.bam ${sample_name}_fixmate.bam ${sample_name}_resorted.bam \\
          ${sample_name}.multimetrics.*.pdf
    """
}

// ============================================================================
// PROCESS: PICARD_METRICS_CRAM  (open-source, metrics-only for pre-aligned CRAM)
// Description: Open-source replacement for SENTIEON_METRICS_CRAM (Ultima path).
//              Runs Picard metrics directly on the pre-aligned CRAM with no
//              dedup re-run (Ultima `demux` already marks duplicates upstream).
//              pct_duplication is derived from the 0x400 flag via samtools flagstat.
//              Picard reads CRAM transparently when given REFERENCE_SEQUENCE.
//              Empty placeholder files match the FASTQ-path glob so the parser
//              needs no separate code path.
//              NOTE: x86-only container (consistent with the FASTQ open-source path).
// ============================================================================
process PICARD_METRICS_CRAM {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(cram), path(crai)
    path fasta_ref
    path base_metrics_intervals
    path wgs_or_target_intervals
    val mode

    output:
    tuple val(sample_name), path("*sentieonmetrics*"), emit: metrics_tuple
    path("*sentieonmetrics*"), emit: metrics_flat

    script:
    def avail_mem = task.memory ? task.memory.giga : 8
    """
    set +u

    # pct_duplication from the existing 0x400 flag (no re-dedup).
    samtools flagstat -@ ${task.cpus} '${cram}' > ${sample_name}.flagstat.sentieonmetrics.txt

    # Picard metric suite — single pass via CollectMultipleMetrics (CRAM read via
    # reference). Runs AlignmentSummary, GcBias, MeanQualityByCycle and
    # QualityYield as PROGRAM modules in ONE pass instead of four. Same Picard
    # collectors → identical values; only file names differ, so we rename to the
    # *sentieonmetrics* names the parser/MultiQC expect. (InsertSize omitted: the
    # CRAM/Ultima path is single-end, matching the original metrics-only set.)
    picard -Xmx${avail_mem}g CollectMultipleMetrics \\
        R=${fasta_ref}/genome.fa \\
        I='${cram}' \\
        O=${sample_name}.multimetrics \\
        PROGRAM=null \\
        PROGRAM=CollectAlignmentSummaryMetrics \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=MeanQualityByCycle \\
        PROGRAM=CollectQualityYieldMetrics

    mv_or_touch() { if [ -f "\$1" ]; then mv "\$1" "\$2"; else touch "\$2"; fi; }
    mv_or_touch ${sample_name}.multimetrics.gc_bias.summary_metrics  ${sample_name}.gcbias_summary.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.gc_bias.detail_metrics   ${sample_name}.gcbias.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.quality_by_cycle_metrics ${sample_name}.meanqualitybycycle.sentieonmetrics.txt
    mv_or_touch ${sample_name}.multimetrics.quality_yield_metrics    ${sample_name}.qualityyield_sentieonmetrics.txt

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
                of.write("\\t".join(parts[:-3]) + "\\n")
            else:
                of.write(line)
PYEOF
    rm -f ${sample_name}.multimetrics.alignment_summary_metrics

    if [ "${mode}" = "exome" ]; then
        # WgsMetrics N/A in exome; hsmetricalgo provided by PICARD_COLLECTHSMETRICS.
        touch ${sample_name}.wgsmetricsalgo.sentieonmetrics.txt
    else
        picard -Xmx${avail_mem}g CollectWgsMetrics \\
            R=${fasta_ref}/genome.fa \\
            I='${cram}' \\
            O=${sample_name}.wgsmetricsalgo.sentieonmetrics.txt \\
            INTERVALS=${wgs_or_target_intervals} \\
            MINIMUM_MAPPING_QUALITY=20 \\
            MINIMUM_BASE_QUALITY=20 \\
            COVERAGE_CAP=250 \\
            COUNT_UNPAIRED=true
        touch ${sample_name}.hsmetricalgo.sentieonmetrics.txt
        touch ${sample_name}.cov_sentieonmetrics
    fi

    rm -f ${sample_name}.multimetrics.*.pdf

    # Placeholders matching the Sentieon CRAM path so the parser glob is uniform.
    # InsertSize is paired-only; dedup is derived from flagstat instead.
    touch ${sample_name}.insertsizemetricalgo.sentieonmetrics.txt
    touch ${sample_name}.dedup_sentieonmetrics.txt
    """
}

// ============================================================================
// PROCESS: WGS_QC_METRICS_TO_PARQUET
// Description: Parse all per-sample QC metrics and write wgsqc_summary Parquet.
//              qc_status is set to PENDING (cluster QC scoring added later).
//              WGS coverage fields are null for exome; HS fields are null for WGS.
//              all_metrics: collected dedup + sentieon + picard (exome only) files.
//              total_reads comes from AlignmentStat (post-subsampling when enabled).
// ============================================================================
process WGS_QC_METRICS_TO_PARQUET {
    tag "wgsqc_metrics_to_parquet"

    input:
    path(all_metrics)
    path(readcount_files)
    val(mode)
    val(mode_prefix)
    val(genome)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    val(user)

    output:
    path("${mode_prefix}qc_summary/workspace=*/workflow_id=*/biosample=*/output.parquet"), emit: parquet
    path("*_${mode_prefix}qc_mqc.txt"), emit: mqc_metrics
    path("${mode_prefix}qc_all_metrics.tsv"), emit: summary_tsv

    script:
    """
    python3 << 'PYEOF'
import os, glob, sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

mode          = "${mode}"
mode_prefix   = "${mode_prefix}"
genome        = "${genome}"
dataset_id    = "${dataset_id}"
workspace     = "${workspace}"
workflow_id   = "${workflow_id}"
pipeline_version = "${pipeline_version}"
user          = "${user}"

# Discover all unique samples from alignment stat files
sample_names = set()
for f in glob.glob("*.alignmentstat_sentieonmetrics.txt"):
    basename = os.path.basename(f)
    sample = basename.split(".alignmentstat")[0]
    sample_names.add(sample)

if not sample_names:
    print("ERROR: No sample metric files found")
    sys.exit(1)

print(f"Found {len(sample_names)} samples: {sorted(sample_names)}")

# ========== Load read counts from CRAM subsampling (Ultima platform only) ==========
# readcount_files are staged by Nextflow into the work directory
readcount_map = {}
for f in glob.glob("*_read_counts.txt"):
    try:
        with open(f, 'r') as fh:
            lines = [line.strip() for line in fh.readlines() if line.strip()]
            if len(lines) >= 1:
                sample = os.path.basename(f).replace("_read_counts.txt", "")
                total = int(lines[0])  # First line is TOTAL_READS
                readcount_map[sample] = {"total_reads": total}
                print(f"Loaded read counts for {sample}: total={total}")
    except Exception as e:
        print(f"Warning: Could not parse {f}: {e}")
print(f"DEBUG: Found {len(readcount_map)} samples with read count files")

all_summaries = []

for sample in sorted(sample_names):
    summary = {
        "biosample":        sample,
        "dataset_id":       dataset_id,
        "pipeline":         "basej-wgs",
        "pipeline_version": pipeline_version,
        "molecule_type":    "dna",
        "mode":             mode,
        "genome":           genome,
        "workspace":        workspace,
        "workflow_id":      workflow_id,
        "user":             user,
        "qc_status":        "PENDING",
    }

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
                    summary["total_reads"]               = int(d.get("TOTAL_READS", 0))
                    summary["pf_reads_aligned"]          = int(d.get("PF_READS_ALIGNED", 0))
                    summary["pf_aligned_bases"]          = int(d.get("PF_ALIGNED_BASES", 0))
                    summary["pf_hq_aligned_reads"]       = int(d.get("PF_HQ_ALIGNED_READS", 0))
                    summary["pf_hq_aligned_bases"]       = int(d.get("PF_HQ_ALIGNED_BASES", 0))
                    summary["pf_hq_aligned_q20_bases"]   = int(d.get("PF_HQ_ALIGNED_Q20_BASES", 0))
                    summary["pf_mismatch_rate"]          = float(d.get("PF_MISMATCH_RATE", 0))
                    summary["pf_hq_error_rate"]          = float(d.get("PF_HQ_ERROR_RATE", 0))
                    summary["pf_indel_rate"]             = float(d.get("PF_INDEL_RATE", 0))
                    summary["mean_read_length"]          = float(d.get("MEAN_READ_LENGTH", 0))
                    summary["pct_reads_aligned_in_pairs"] = float(d.get("PCT_READS_ALIGNED_IN_PAIRS", 0))
                    summary["pct_chimeras"]              = float(d.get("PCT_CHIMERAS", 0)) if d.get("PCT_CHIMERAS") else 0.0
                    summary["pct_adapter"]               = float(d.get("PCT_ADAPTER", 0))
                    summary["strand_balance"]            = float(d.get("STRAND_BALANCE", 0))
        except Exception as e:
            print(f"Warning: Could not parse alignment stats for {sample}: {e}")

    # ========== Parse Sentieon QualityYield ==========
    # Emits Q20_BASES / Q30_BASES / TOTAL_BASES so we can report a true Q30 rate
    # (matches what fastp would compute for FASTQ input). Falls back gracefully
    # if the file is absent.
    for f in glob.glob(f"{sample}*.qualityyield_sentieonmetrics.txt"):
        if os.path.getsize(f) == 0:
            continue
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
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
                        break
        except Exception as e:
            print(f"Warning: Could not parse QualityYield for {sample}: {e}")

    # ========== Override total_reads from CRAM subsampling if available ==========
    if sample in readcount_map:
        summary["total_reads"] = readcount_map[sample]["total_reads"]
        print(f"Overrode total_reads for {sample} with CRAM subsampling count: {summary['total_reads']}")

    # ========== Parse Dedup Metrics ==========
    for f in glob.glob(f"{sample}*.dedup_sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                for i, line in enumerate(lines):
                    if line.startswith("LIBRARY"):
                        headers = line.strip().split("\\t")
                        if i + 1 < len(lines):
                            parts = lines[i + 1].strip().split("\\t")
                            d = dict(zip(headers, parts))
                            summary["pct_duplication"]            = float(d.get("PERCENT_DUPLICATION", 0))
                            summary["estimated_library_size"]     = int(d.get("ESTIMATED_LIBRARY_SIZE", 0)) if d.get("ESTIMATED_LIBRARY_SIZE") else 0
                            summary["read_pairs_examined"]        = int(d.get("READ_PAIRS_EXAMINED", 0))
                            summary["read_pair_duplicates"]       = int(d.get("READ_PAIR_DUPLICATES", 0))
                            summary["read_pair_optical_duplicates"] = int(d.get("READ_PAIR_OPTICAL_DUPLICATES", 0))
                        break
        except Exception as e:
            print(f"Warning: Could not parse dedup for {sample}: {e}")

    # ========== Fallback: pct_duplication from samtools flagstat (CRAM path) ==========
    # When Sentieon Dedup did not run (CRAM input where duplicates are already
    # marked upstream, e.g. Ultima `demux`), derive PERCENT_DUPLICATION from
    # the existing 0x400 flag via samtools flagstat. Only applied when the
    # FASTQ-path dedup parser above did not already populate the field.
    #
    # Sample flagstat lines (samtools >=1.13):
    #   12000 + 0 primary
    #   500   + 0 primary duplicates
    #   10500 + 0 primary mapped (95.45% : N/A)
    # We use the QC-pass count (first integer on each line).
    if "pct_duplication" not in summary:
        for f in glob.glob(f"{sample}*.flagstat.sentieonmetrics.txt"):
            if os.path.getsize(f) == 0:
                continue
            try:
                primary = primary_dup = primary_mapped = 0
                with open(f, 'r') as fh:
                    for line in fh:
                        s = line.strip()
                        if not s:
                            continue
                        toks = s.split()
                        try:
                            n = int(toks[0])
                        except (ValueError, IndexError):
                            continue
                        # Strip the leading "N + M " count prefix
                        # to leave just the descriptor (e.g. "primary",
                        # "primary duplicates", "primary mapped (...)").
                        rest = s.split(None, 3)[-1] if len(toks) > 3 else ""
                        if rest.startswith("primary mapped"):
                            primary_mapped = n
                        elif rest.startswith("primary duplicates"):
                            primary_dup = n
                        elif rest == "primary":
                            primary = n
                # Prefer mapped denominator (matches Picard PERCENT_DUPLICATION
                # which counts duplicates among mapped reads). Fall back to
                # total primary if mapped is unavailable.
                denom = primary_mapped if primary_mapped > 0 else primary
                summary["pct_duplication"] = (primary_dup / denom) if denom > 0 else 0.0
                summary["read_pairs_examined"] = denom
                summary["read_pair_duplicates"] = primary_dup
            except Exception as e:
                print(f"Warning: Could not parse flagstat for {sample}: {e}")

    # ========== Parse Insert Size Metrics ==========
    for f in glob.glob(f"{sample}*.insertsizemetricalgo.sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                for i, line in enumerate(lines):
                    if line.startswith("MEDIAN_INSERT_SIZE"):
                        headers = line.strip().split("\\t")
                        if i + 1 < len(lines):
                            parts = lines[i + 1].strip().split("\\t")
                            d = dict(zip(headers, parts))
                            summary["insert_median"]    = float(d.get("MEDIAN_INSERT_SIZE", 0))
                            summary["insert_mad"]       = float(d.get("MEDIAN_ABSOLUTE_DEVIATION", 0))
                            summary["insert_min"]       = int(d.get("MIN_INSERT_SIZE", 0)) if d.get("MIN_INSERT_SIZE") else 0
                            summary["insert_max"]       = int(d.get("MAX_INSERT_SIZE", 0)) if d.get("MAX_INSERT_SIZE") else 0
                            summary["insert_mean"]      = float(d.get("MEAN_INSERT_SIZE", 0))
                            summary["insert_std"]       = float(d.get("STANDARD_DEVIATION", 0))
                            summary["insert_read_pairs"] = int(d.get("READ_PAIRS", 0)) if d.get("READ_PAIRS") else 0
                        break
        except Exception as e:
            print(f"Warning: Could not parse insert size for {sample}: {e}")

    # ========== Parse GC Bias Summary ==========
    for f in glob.glob(f"{sample}*.gcbias_summary.sentieonmetrics.txt"):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
                for i, line in enumerate(lines):
                    if line.startswith("ACCUMULATION"):
                        headers = line.strip().split("\\t")
                        if i + 1 < len(lines):
                            parts = lines[i + 1].strip().split("\\t")
                            d = dict(zip(headers, parts))
                            summary["at_dropout"]    = float(d.get("AT_DROPOUT", 0))
                            summary["gc_dropout"]    = float(d.get("GC_DROPOUT", 0))
                            summary["gc_nc_0_19"]    = float(d.get("GC_NC_0_19", 0))
                            summary["gc_nc_20_39"]   = float(d.get("GC_NC_20_39", 0))
                            summary["gc_nc_40_59"]   = float(d.get("GC_NC_40_59", 0))
                            summary["gc_nc_60_79"]   = float(d.get("GC_NC_60_79", 0))
                            summary["gc_nc_80_100"]  = float(d.get("GC_NC_80_100", 0))
                        break
        except Exception as e:
            print(f"Warning: Could not parse GC bias for {sample}: {e}")

    # ========== Initialize WGS / HS fields as None ==========
    summary["mean_coverage"]    = None
    summary["sd_coverage"]      = None
    summary["median_coverage"]  = None
    summary["mad_coverage"]     = None
    summary["pct_1x"]           = None
    summary["pct_5x"]           = None
    summary["pct_10x"]          = None
    summary["pct_15x"]          = None
    summary["pct_20x"]          = None
    summary["pct_25x"]          = None
    summary["pct_30x"]          = None
    summary["pct_40x"]          = None
    summary["pct_50x"]          = None
    summary["pct_60x"]          = None
    summary["pct_70x"]          = None
    summary["pct_80x"]          = None
    summary["pct_90x"]          = None
    summary["pct_100x"]         = None
    summary["pct_exc_mapq"]     = None
    summary["pct_exc_dupe"]     = None
    summary["pct_exc_unpaired"] = None
    summary["pct_exc_baseq"]    = None
    summary["pct_exc_overlap"]  = None
    summary["pct_exc_capped"]   = None
    summary["pct_exc_total"]    = None
    summary["het_snp_sensitivity"] = None
    summary["het_snp_q"]        = None

    summary["mean_target_coverage"]  = None
    summary["fold_enrichment"]       = None
    summary["pct_selected_bases"]    = None
    summary["fold_80_base_penalty"]  = None
    summary["zero_cvg_targets_pct"]  = None
    summary["on_bait_bases"]         = None
    summary["near_bait_bases"]       = None
    summary["off_bait_bases"]        = None
    summary["pct_target_bases_1x"]    = None
    summary["pct_target_bases_2x"]   = None
    summary["pct_target_bases_10x"]  = None
    summary["pct_target_bases_20x"]  = None
    summary["pct_target_bases_30x"]  = None
    summary["pct_target_bases_40x"]  = None
    summary["pct_target_bases_50x"]  = None
    summary["pct_target_bases_100x"] = None
    summary["pct_target_bases_250x"] = None
    summary["pct_target_bases_500x"] = None
    summary["pct_target_bases_1000x"] = None

    # ========== Parse WGS Metrics (WGS mode only) ==========
    if mode == 'wgs':
        for f in glob.glob(f"{sample}*.wgsmetricsalgo.sentieonmetrics.txt"):
            if os.path.getsize(f) == 0:
                continue
            try:
                with open(f, 'r') as fh:
                    lines = fh.readlines()
                    for i, line in enumerate(lines):
                        if line.startswith("GENOME_TERRITORY") or line.startswith("MEAN_COVERAGE"):
                            headers = line.strip().split("\\t")
                            if i + 1 < len(lines):
                                parts = lines[i + 1].strip().split("\\t")
                                d = dict(zip(headers, parts))
                                summary["mean_coverage"]    = float(d.get("MEAN_COVERAGE", 0))
                                summary["sd_coverage"]      = float(d.get("SD_COVERAGE", 0))
                                summary["median_coverage"]  = float(d.get("MEDIAN_COVERAGE", 0))
                                summary["mad_coverage"]     = float(d.get("MAD_COVERAGE", 0))
                                summary["pct_1x"]           = float(d.get("PCT_1X", 0))
                                summary["pct_5x"]           = float(d.get("PCT_5X", 0))
                                summary["pct_10x"]          = float(d.get("PCT_10X", 0))
                                summary["pct_15x"]          = float(d.get("PCT_15X", 0))
                                summary["pct_20x"]          = float(d.get("PCT_20X", 0))
                                summary["pct_25x"]          = float(d.get("PCT_25X", 0))
                                summary["pct_30x"]          = float(d.get("PCT_30X", 0))
                                summary["pct_40x"]          = float(d.get("PCT_40X", 0))
                                summary["pct_50x"]          = float(d.get("PCT_50X", 0))
                                summary["pct_60x"]          = float(d.get("PCT_60X", 0))
                                summary["pct_70x"]          = float(d.get("PCT_70X", 0))
                                summary["pct_80x"]          = float(d.get("PCT_80X", 0))
                                summary["pct_90x"]          = float(d.get("PCT_90X", 0))
                                summary["pct_100x"]         = float(d.get("PCT_100X", 0))
                                summary["pct_exc_mapq"]     = float(d.get("PCT_EXC_MAPQ", 0))
                                summary["pct_exc_dupe"]     = float(d.get("PCT_EXC_DUPE", 0))
                                summary["pct_exc_unpaired"] = float(d.get("PCT_EXC_UNPAIRED", 0))
                                summary["pct_exc_baseq"]    = float(d.get("PCT_EXC_BASEQ", 0))
                                summary["pct_exc_overlap"]  = float(d.get("PCT_EXC_OVERLAP", 0))
                                summary["pct_exc_capped"]   = float(d.get("PCT_EXC_CAPPED", 0))
                                summary["pct_exc_total"]    = float(d.get("PCT_EXC_TOTAL", 0))
                                summary["het_snp_sensitivity"] = float(d.get("HET_SNP_SENSITIVITY", 0))
                                summary["het_snp_q"]        = int(d.get("HET_SNP_Q", 0)) if d.get("HET_SNP_Q") else 0
                            break
            except Exception as e:
                print(f"Warning: Could not parse WGS metrics for {sample}: {e}")

    # ========== Parse HS Metrics (exome mode only) ==========
    if mode == 'exome':
        for f in glob.glob(f"{sample}*.hsmetricalgo.sentieonmetrics.txt"):
            if os.path.getsize(f) == 0:
                continue
            try:
                with open(f, 'r') as fh:
                    lines = fh.readlines()
                    for i, line in enumerate(lines):
                        if "MEAN_TARGET_COVERAGE" in line and not line.startswith("#"):
                            headers = line.strip().split("\\t")
                            if i + 1 < len(lines):
                                parts = lines[i + 1].strip().split("\\t")
                                d = dict(zip(headers, parts))
                                summary["mean_target_coverage"]  = float(d.get("MEAN_TARGET_COVERAGE", 0))
                                summary["fold_enrichment"]       = float(d.get("FOLD_ENRICHMENT", 0))
                                summary["pct_selected_bases"]    = float(d.get("PCT_SELECTED_BASES", 0))
                                raw_f80 = d.get("FOLD_80_BASE_PENALTY", None)
                                summary["fold_80_base_penalty"]  = float(raw_f80) if raw_f80 and raw_f80 != '?' else None
                                summary["zero_cvg_targets_pct"]  = float(d.get("ZERO_CVG_TARGETS_PCT", 0))
                                summary["on_bait_bases"]         = int(d.get("ON_BAIT_BASES", 0))
                                summary["near_bait_bases"]       = int(d.get("NEAR_BAIT_BASES", 0))
                                summary["off_bait_bases"]        = int(d.get("OFF_BAIT_BASES", 0))
                                summary["pct_target_bases_1x"]    = float(d.get("PCT_TARGET_BASES_1X", 0))
                                summary["pct_target_bases_2x"]   = float(d.get("PCT_TARGET_BASES_2X", 0))
                                summary["pct_target_bases_10x"]  = float(d.get("PCT_TARGET_BASES_10X", 0))
                                summary["pct_target_bases_20x"]  = float(d.get("PCT_TARGET_BASES_20X", 0))
                                summary["pct_target_bases_30x"]  = float(d.get("PCT_TARGET_BASES_30X", 0))
                                summary["pct_target_bases_40x"]  = float(d.get("PCT_TARGET_BASES_40X", 0))
                                summary["pct_target_bases_50x"]  = float(d.get("PCT_TARGET_BASES_50X", 0))
                                summary["pct_target_bases_100x"] = float(d.get("PCT_TARGET_BASES_100X", 0))
                                summary["pct_target_bases_250x"] = float(d.get("PCT_TARGET_BASES_250X", 0))
                                summary["pct_target_bases_500x"] = float(d.get("PCT_TARGET_BASES_500X", 0))
                                summary["pct_target_bases_1000x"] = float(d.get("PCT_TARGET_BASES_1000X", 0))
                            break
            except Exception as e:
                print(f"Warning: Could not parse HS metrics for {sample}: {e}")

    # ========== Ensure all expected fields are set ==========
    if "pct_chimeras" not in summary:
        summary["pct_chimeras"] = 0.0
    if "total_reads" not in summary:
        summary["total_reads"] = 0

    all_summaries.append(summary)

    # ========== Write MultiQC TSV ==========
    with open(f"{sample}_{mode_prefix}qc_mqc.txt", 'w') as f:
        section_desc = 'Per-sample alignment, coverage, and duplication QC metrics with pass/fail status.' if mode == 'wgs' else 'Per-sample alignment, target coverage, and duplication QC metrics with pass/fail status.'
        f.write(f"# id: '{mode_prefix}qc_summary'\\n")
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'QC Summary'\\n")
        f.write(f"# description: '{section_desc}'\\n")
        f.write("# pconfig:\\n")
        f.write(f"#   id: '{mode_prefix}qc_summary_table'\\n")
        if mode == 'wgs':
            f.write("sample_name\\tQC_Status\\tScore\\tTotal_Reads\\tPCT_Duplication\\tPCT_1x\\tPCT_5x\\tPCT_10x\\tPCT_30x\\tPCT_Chimeras\\tMean_Coverage\\tInsert_Median\\tAT_Dropout\\tGC_Dropout\\n")
            f.write(f"{sample}\\tPENDING\\tNA\\t{summary.get('total_reads', 0)}\\t{summary.get('pct_duplication', 0):.4f}\\t")
            f.write(f"{summary.get('pct_1x', 0) or 0:.4f}\\t{summary.get('pct_5x', 0) or 0:.4f}\\t")
            f.write(f"{summary.get('pct_10x', 0) or 0:.4f}\\t{summary.get('pct_30x', 0) or 0:.4f}\\t")
            f.write(f"{summary.get('pct_chimeras', 0):.4f}\\t{summary.get('mean_coverage', 0) or 0:.2f}\\t")
            f.write(f"{summary.get('insert_median', 0) or 0:.0f}\\t{summary.get('at_dropout', 0) or 0:.4f}\\t{summary.get('gc_dropout', 0) or 0:.4f}\\n")
        else:
            f.write("sample_name\\tQC_Status\\tScore\\tTotal_Reads\\tPCT_Target_10x\\tZero_Cvg_Targets_Pct\\tFold_80_Base_Penalty\\tMean_Target_Coverage\\tPCT_Selected_Bases\\tInsert_Median\\tPCT_Chimeras\\n")
            f.write(f"{sample}\\tPENDING\\tNA\\t{summary.get('total_reads', 0)}\\t")
            f.write(f"{summary.get('pct_target_bases_10x', 0) or 0:.4f}\\t{summary.get('zero_cvg_targets_pct', 0) or 0:.4f}\\t")
            f.write(f"{summary.get('fold_80_base_penalty', 0) or 0:.4f}\\t{summary.get('mean_target_coverage', 0) or 0:.2f}\\t")
            f.write(f"{summary.get('pct_selected_bases', 0) or 0:.4f}\\t{summary.get('insert_median', 0) or 0:.0f}\\t")
            f.write(f"{summary.get('pct_chimeras', 0):.4f}\\n")

    print(f"  {sample}: metrics parsed, qc_status=PENDING")

# ========== Write Parquets ==========
_str_cols    = ['biosample','dataset_id','pipeline','pipeline_version','molecule_type',
                'mode','genome','workspace','workflow_id','user','qc_status']
_double_cols = ['pct_duplication','pf_mismatch_rate','pf_hq_error_rate','pf_indel_rate',
                'pct_reads_aligned_in_pairs','pct_chimeras','pct_adapter','strand_balance',
                'mean_read_length',
                'sentieon_q20_rate','sentieon_q30_rate',
                'insert_median','insert_mad','insert_mean','insert_std',
                'at_dropout','gc_dropout',
                'gc_nc_0_19','gc_nc_20_39','gc_nc_40_59','gc_nc_60_79','gc_nc_80_100',
                'mean_coverage','sd_coverage','median_coverage','mad_coverage',
                'pct_1x','pct_5x','pct_10x','pct_15x','pct_20x','pct_25x','pct_30x',
                'pct_40x','pct_50x','pct_60x','pct_70x','pct_80x','pct_90x','pct_100x',
                'pct_exc_mapq','pct_exc_dupe','pct_exc_unpaired','pct_exc_baseq',
                'pct_exc_overlap','pct_exc_capped','pct_exc_total',
                'het_snp_sensitivity',
                'mean_target_coverage','fold_enrichment','pct_selected_bases',
                'fold_80_base_penalty','zero_cvg_targets_pct',
                'pct_target_bases_1x','pct_target_bases_2x','pct_target_bases_10x','pct_target_bases_20x',
                'pct_target_bases_30x','pct_target_bases_40x','pct_target_bases_50x','pct_target_bases_100x',
                'pct_target_bases_250x','pct_target_bases_500x','pct_target_bases_1000x']
_bigint_cols = ['total_reads','pf_reads_aligned',
                'pf_aligned_bases','pf_hq_aligned_reads','pf_hq_aligned_bases','pf_hq_aligned_q20_bases',
                'estimated_library_size','read_pairs_examined',
                'read_pair_duplicates','read_pair_optical_duplicates','insert_read_pairs',
                'insert_min','insert_max',
                'sentieon_total_bases',
                'het_snp_q',
                'on_bait_bases','near_bait_bases','off_bait_bases']

for summary in all_summaries:
    sample = summary["biosample"]
    df = pd.DataFrame([summary])
    for col in _double_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').astype('float64')
    for col in _bigint_cols:
        if col in df.columns:
            s = pd.to_numeric(df[col], errors='coerce')
            df[col] = s.where(s.isna(), s.round(0)).astype('Int64')
    for col in _str_cols:
        if col in df.columns:
            df[col] = df[col].astype('string')
    out_dir = f"{mode_prefix}qc_summary/workspace={workspace}/workflow_id={workflow_id}/biosample={sample}"
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   os.path.join(out_dir, "output.parquet"))
    print(f"  {sample}: parquet written")

# ========== Write combined TSV ==========
pd.DataFrame(all_summaries).to_csv(f"{mode_prefix}qc_all_metrics.tsv", sep="\\t", index=False)
print(f"Written {len(all_summaries)} parquets and combined TSV")
PYEOF
    """
}

// ============================================================================
// PROCESS: WGS_QC_PLOTS
// Description: Run mode-aware R QC clustering/scoring script (wgs_qc_plot.R
//              for wgs, wes_qc_plot.R for exome), then rewrite parquets and
//              all_metrics TSV with final qc_status and qc_score.
//              Thresholds — WGS (max 5): >=4 PASS, 3 Borderline, <3 FAIL.
//                           WES (max 4): >=3 PASS, 2 Borderline, <2 FAIL.
// ============================================================================
process WGS_QC_PLOTS {
    tag "wgs_qc_plots"

    input:
    path(all_metrics_tsv)
    val(mode)
    val(mode_prefix)
    val(workspace)
    val(workflow_id)

    output:
    path("${mode_prefix}qc_summary/workspace=*/workflow_id=*/biosample=*/output.parquet"), emit: parquet
    path("${mode_prefix}qc_all_metrics.tsv"),                                              emit: summary_tsv
    path("*.pdf"),                                                                          emit: plots_pdf
    path("*_mqc.jpg"),                                                                      emit: plots_jpg
    path("*_${mode_prefix}qc_mqc.txt"),                                                    emit: mqc_metrics
    path("*-QC_ConsensusScores.txt"),                                                       emit: scores
    path("*-QC_ConsensusScores_SummaryTable_mqc.txt"),                                      emit: scores_summary
    path("per_biosample_status.csv"),                                                       emit: per_biosample_status

    script:
    def r_script    = mode == 'exome' ? 'wes_qc_plot.R' : 'wgs_qc_plot.R'
    def scores_file = mode == 'exome' ? 'WES-QC_ConsensusScores.txt' : 'WGS-QC_ConsensusScores.txt'
    """
    # Make a physical copy so we can write a fresh ${mode_prefix}qc_all_metrics.tsv
    # without following the staged symlink back to the previous process work dir
    cp ${all_metrics_tsv} input_metrics.tsv
    rm -f ${mode_prefix}qc_all_metrics.tsv

    # Run mode-aware R QC plot script
    Rscript /usr/local/bin/${r_script} --metrics_file input_metrics.tsv

    # Update parquets and all_metrics TSV with final QC status from R scores
    python3 << 'PYEOF'
import os, sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

mode        = "${mode}"
mode_prefix = "${mode_prefix}"
workspace   = "${workspace}"
workflow_id = "${workflow_id}"
scores_file = "${scores_file}"

# ========== Read R consensus scores ==========
if not os.path.exists(scores_file):
    print(f"WARNING: {scores_file} not found; all samples will keep qc_status=PENDING")
    scores_df = pd.DataFrame(columns=["biosample", "qc_score", "qc_status"])
else:
    scores_df = pd.read_csv(scores_file, sep="\\t")
    scores_df = scores_df[["SampleId", "CompositeScore"]].rename(
        columns={"SampleId": "biosample", "CompositeScore": "qc_score"})
    scores_df["qc_score"] = pd.to_numeric(scores_df["qc_score"], errors="coerce")

    def score_to_status(s):
        if pd.isna(s):
            return "PENDING"
        s = int(s)
        if mode == "wgs":
            if s >= 4: return "PASS"
            if s == 3: return "Borderline"
            return "FAIL"
        else:  # exome
            if s >= 3: return "PASS"
            if s == 2: return "Borderline"
            return "FAIL"

    scores_df["qc_status"] = scores_df["qc_score"].apply(score_to_status)

score_map  = dict(zip(scores_df["biosample"], scores_df["qc_score"]))
status_map = dict(zip(scores_df["biosample"], scores_df["qc_status"]))

# ========== Read pending all_metrics TSV ==========
all_df = pd.read_csv("input_metrics.tsv", sep="\\t")

# Update qc_status and add qc_score column
all_df["qc_status"] = all_df["biosample"].map(status_map).fillna("PENDING")
all_df["qc_score"]  = all_df["biosample"].map(score_map)

# Override score for zero-read samples — R scoring assigns non-zero scores
# to samples with all-zero metrics, but they should show as FAIL / NA
for idx, row in all_df.iterrows():
    if int(row.get("total_reads", 0) or 0) == 0:
        all_df.at[idx, "qc_status"] = "FAIL"
        all_df.at[idx, "qc_score"]  = None
        score_map[row["biosample"]]  = None
        status_map[row["biosample"]] = "FAIL"

print(f"Found {len(all_df)} samples to update")

# ========== Schema definitions (mirrors WGS_QC_METRICS_TO_PARQUET + qc_score) ==========
_str_cols    = ['biosample', 'dataset_id', 'pipeline', 'pipeline_version', 'molecule_type',
                'mode', 'genome', 'workspace', 'workflow_id', 'user', 'qc_status']
_double_cols = ['pct_duplication', 'pf_mismatch_rate', 'pf_hq_error_rate', 'pf_indel_rate',
                'pct_reads_aligned_in_pairs', 'pct_chimeras', 'pct_adapter', 'strand_balance',
                'mean_read_length',
                'sentieon_q20_rate', 'sentieon_q30_rate',
                'insert_median', 'insert_mad', 'insert_mean', 'insert_std',
                'at_dropout', 'gc_dropout',
                'gc_nc_0_19', 'gc_nc_20_39', 'gc_nc_40_59', 'gc_nc_60_79', 'gc_nc_80_100',
                'mean_coverage', 'sd_coverage', 'median_coverage', 'mad_coverage',
                'pct_1x', 'pct_5x', 'pct_10x', 'pct_15x', 'pct_20x', 'pct_25x', 'pct_30x',
                'pct_40x', 'pct_50x', 'pct_60x', 'pct_70x', 'pct_80x', 'pct_90x', 'pct_100x',
                'pct_exc_mapq', 'pct_exc_dupe', 'pct_exc_unpaired', 'pct_exc_baseq',
                'pct_exc_overlap', 'pct_exc_capped', 'pct_exc_total',
                'het_snp_sensitivity',
                'mean_target_coverage', 'fold_enrichment', 'pct_selected_bases',
                'fold_80_base_penalty', 'zero_cvg_targets_pct',
                'pct_target_bases_1x', 'pct_target_bases_2x', 'pct_target_bases_10x', 'pct_target_bases_20x',
                'pct_target_bases_30x', 'pct_target_bases_40x', 'pct_target_bases_50x', 'pct_target_bases_100x',
                'pct_target_bases_250x', 'pct_target_bases_500x', 'pct_target_bases_1000x']
_bigint_cols = ['total_reads', 'pf_reads_aligned',
                'pf_aligned_bases', 'pf_hq_aligned_reads', 'pf_hq_aligned_bases', 'pf_hq_aligned_q20_bases',
                'estimated_library_size', 'read_pairs_examined',
                'read_pair_duplicates', 'read_pair_optical_duplicates', 'insert_read_pairs',
                'insert_min', 'insert_max',
                'sentieon_total_bases',
                'het_snp_q',
                'on_bait_bases', 'near_bait_bases', 'off_bait_bases',
                'qc_score']

# ========== Write per-sample final parquets ==========
for _, row in all_df.iterrows():
    sample = row["biosample"]
    df = pd.DataFrame([row.to_dict()])
    for col in _double_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').astype('float64')
    for col in _bigint_cols:
        if col in df.columns:
            s = pd.to_numeric(df[col], errors='coerce')
            df[col] = s.where(s.isna(), s.round(0)).astype('Int64')
    for col in _str_cols:
        if col in df.columns:
            df[col] = df[col].astype('string')
    out_dir = f"{mode_prefix}qc_summary/workspace={workspace}/workflow_id={workflow_id}/biosample={sample}"
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   os.path.join(out_dir, "output.parquet"))
    print(f"  {sample}: qc_status={status_map.get(sample, 'PENDING')}, qc_score={score_map.get(sample)}")

# ========== Write updated all_metrics TSV ==========
all_df.to_csv(f"{mode_prefix}qc_all_metrics.tsv", sep="\\t", index=False)
all_df[['biosample','qc_status','pipeline','pipeline_version']].rename(
    columns={'biosample':'biosampleName'}).to_csv('per_biosample_status.csv', index=False)
print(f"Written {len(all_df)} final parquets and updated all_metrics TSV")

# ========== Rewrite per-sample mqc files with final QC status and Score ==========
for _, row in all_df.iterrows():
    sample = row["biosample"]
    status = row.get("qc_status", "PENDING")
    score_val = score_map.get(sample)
    score_str = f"{int(score_val)}" if score_val is not None and not pd.isna(score_val) else "NA"
    with open(f"{sample}_{mode_prefix}qc_mqc.txt", 'w') as f:
        section_desc = 'Per-sample alignment, coverage, and duplication QC metrics with pass/fail status.' if mode == 'wgs' else 'Per-sample alignment, target coverage, and duplication QC metrics with pass/fail status.'
        f.write(f"# id: '{mode_prefix}qc_summary'\\n")
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'QC Summary'\\n")
        f.write(f"# description: '{section_desc}'\\n")
        f.write("# pconfig:\\n")
        f.write(f"#   id: '{mode_prefix}qc_summary_table'\\n")
        if mode == "wgs":
            f.write("sample_name\\tQC_Status\\tScore\\tTotal_Reads\\tPCT_Duplication\\tPCT_1x\\tPCT_5x\\tPCT_10x\\tPCT_30x\\tPCT_Chimeras\\tMean_Coverage\\tInsert_Median\\tAT_Dropout\\tGC_Dropout\\n")
            f.write(f"{sample}\\t{status}\\t{score_str}\\t{int(row.get('total_reads', 0) or 0)}\\t{float(row.get('pct_duplication', 0) or 0):.4f}\\t")
            f.write(f"{float(row.get('pct_1x', 0) or 0):.4f}\\t{float(row.get('pct_5x', 0) or 0):.4f}\\t")
            f.write(f"{float(row.get('pct_10x', 0) or 0):.4f}\\t{float(row.get('pct_30x', 0) or 0):.4f}\\t")
            f.write(f"{float(row.get('pct_chimeras', 0) or 0):.4f}\\t{float(row.get('mean_coverage', 0) or 0):.2f}\\t")
            f.write(f"{float(row.get('insert_median', 0) or 0):.0f}\\t{float(row.get('at_dropout', 0) or 0):.4f}\\t{float(row.get('gc_dropout', 0) or 0):.4f}\\n")
        else:
            f.write("sample_name\\tQC_Status\\tScore\\tTotal_Reads\\tPCT_Target_10x\\tZero_Cvg_Targets_Pct\\tFold_80_Base_Penalty\\tMean_Target_Coverage\\tPCT_Selected_Bases\\tInsert_Median\\tPCT_Chimeras\\n")
            f.write(f"{sample}\\t{status}\\t{score_str}\\t{int(row.get('total_reads', 0) or 0)}\\t")
            f.write(f"{float(row.get('pct_target_bases_10x', 0) or 0):.4f}\\t{float(row.get('zero_cvg_targets_pct', 0) or 0):.4f}\\t")
            f.write(f"{float(row.get('fold_80_base_penalty', 0) or 0):.4f}\\t{float(row.get('mean_target_coverage', 0) or 0):.2f}\\t")
            f.write(f"{float(row.get('pct_selected_bases', 0) or 0):.4f}\\t{float(row.get('insert_median', 0) or 0):.0f}\\t")
            f.write(f"{float(row.get('pct_chimeras', 0) or 0):.4f}\\n")
print(f"Rewrote {len(all_df)} mqc files with final QC status and Score")
PYEOF

    # Prepend MultiQC custom-content headers so the SummaryTable renders as a table in the report
    summary_table=\$(ls *-QC_ConsensusScores_SummaryTable_mqc.txt)
    tmp=\$(mktemp)
    printf '# id: "qc_score_distribution"\n# section_name: "Total Usable Cells"\n# description: "Distribution of samples across composite QC score categories. Scores range from 0 (poor quality) to 5 (high quality)."\n# plot_type: "table"\n# pconfig:\n#   id: "qc_score_dist_table"\n#   title: "Total Usable Cells"\n' | cat - "\$summary_table" > "\$tmp" && mv "\$tmp" "\$summary_table"
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
    val(mode)
    val(dataset_id)
    val(workspace)
    val(workflow_id)
    val(pipeline_version)
    path(logo)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    def report_title = mode == 'exome' ? "basej-wesqc v${pipeline_version}" : "basej-wgsqc v${pipeline_version}"
    def mode_prefix  = mode == 'exome' ? 'wes'                              : 'wgs'
    """
    cat > multiqc_config.yaml << EOF
custom_logo_title: 'BioSkryb Genomics'
custom_logo: bioskryb_logo-tagline.png
custom_logo_width: 260

title: "${report_title}"
report_header_info:
  - Dataset ID: "${dataset_id}"
  - Workspace: "${workspace}"
  - Workflow ID: "${workflow_id}"
show_analysis_paths: false
show_analysis_time: false
skip_generalstats: true

fn_clean_exts:
  - ".fastq.gz"
  - ".fq.gz"
  - ".bam"
  - "_sorted"
  - "_sentieonmetrics"
  - ".txt"
  - ".json"

extra_fn_clean_exts:
  - "_R1"
  - "_R2"
  - "_sorted"
  - "_sentieonmetrics"
  - ".dedup"

module_order:
  - custom_content
  - picard

report_section_order:
  ${mode_prefix}qc_summary:
    order: 1000
  ${mode_prefix.toUpperCase()}-QC_composition_mqc.jpg:
    order: 900
  qc_score_distribution:
    order: 800

table_columns_visible:
  picard:
    PERCENT_DUPLICATION: true
    summed_median: true
    MEAN_TARGET_COVERAGE: true
    PCT_TARGET_BASES_20X: true

table_cond_formatting_rules:
  QC_Status:
    pass:
      - s_eq: "PASS"
    warn:
      - s_eq: "Borderline"
    fail:
      - s_eq: "FAIL"

custom_data:
  ${mode_prefix.toUpperCase()}-QC_composition_mqc.jpg:
    section_name: "QC Composition"
    description: "Per-sample QC metric distributions grouped by cluster. Each panel shows a key metric across all cells."
  ${mode_prefix}qc_summary:
    headers:
      QC_Status:
        title: "QC Status"
        description: "Composite QC status: PASS, Borderline, or FAIL"
        placement: 100
      Score:
        title: "Score"
        description: "Composite QC score as an integer (max 5 for WGS, 4 for exome)"
        placement: 110
      Total_Reads:
        title: "Total Reads"
        description: "Total number of read pairs"
        format: "{:,.0f}"
        placement: 120
      PCT_Duplication:
        title: "Proportion Duplication"
        description: "Proportion of duplicate read pairs (0-1). Source: Sentieon Dedup PERCENT_DUPLICATION. Note: coordinate-based detection may undercount duplicates in amplification-based chemistries."
        format: "{:.4f}"
        placement: 130
      PCT_1x:
        title: "Proportion >=1x"
        description: "Proportion of genome bases covered at >=1x depth (0-1). Source: Sentieon WgsMetricsAlgo PCT_1X."
        format: "{:.4f}"
        placement: 140
      PCT_5x:
        title: "Proportion >=5x"
        description: "Proportion of genome bases covered at >=5x depth (0-1). Source: Sentieon WgsMetricsAlgo PCT_5X."
        format: "{:.4f}"
        placement: 150
      PCT_10x:
        title: "Proportion >=10x"
        description: "Proportion of genome bases covered at >=10x depth (0-1). Source: Sentieon WgsMetricsAlgo PCT_10X."
        format: "{:.4f}"
        placement: 160
      PCT_30x:
        title: "Proportion >=30x"
        description: "Proportion of genome bases covered at >=30x depth (0-1). Source: Sentieon WgsMetricsAlgo PCT_30X."
        format: "{:.4f}"
        placement: 170
      PCT_Chimeras:
        title: "Proportion Chimeras"
        description: "Proportion of chimeric read pairs (0-1). Source: Sentieon AlignmentStat PCT_CHIMERAS."
        format: "{:.4f}"
        placement: 180
      Mean_Coverage:
        title: "Mean Coverage"
        description: "Mean genome coverage depth (x). Source: Sentieon WgsMetricsAlgo MEAN_COVERAGE."
        format: "{:.2f}"
        placement: 190
      Insert_Median:
        title: "Insert Median (bp)"
        description: "Median insert size in base pairs. Source: Sentieon InsertSizeMetricAlgo."
        format: "{:.0f}"
        placement: 200
      AT_Dropout:
        title: "AT Dropout"
        description: "AT dropout score — measure of coverage loss in AT-rich regions. Source: Sentieon WgsMetricsAlgo AT_DROPOUT."
        format: "{:.4f}"
        placement: 210
      GC_Dropout:
        title: "GC Dropout"
        description: "GC dropout score — measure of coverage loss in GC-rich regions. Source: Sentieon WgsMetricsAlgo GC_DROPOUT."
        format: "{:.4f}"
        placement: 212
      PCT_Target_10x:
        title: "Proportion Target >=10x"
        description: "Proportion of target bases covered at >=10x depth (0-1). Source: Picard CollectHsMetrics PCT_TARGET_BASES_10X."
        format: "{:.4f}"
        placement: 140
      Zero_Cvg_Targets_Pct:
        title: "Proportion Zero Coverage Targets"
        description: "Proportion of target bases with zero coverage (0-1). Source: Picard CollectHsMetrics ZERO_CVG_TARGETS_PCT."
        format: "{:.4f}"
        placement: 150
      Fold_80_Base_Penalty:
        title: "Fold 80 Base Penalty"
        description: "Fold 80 base penalty — ratio of coverage needed to cover 80% of bases vs mean coverage. Lower = more uniform. Source: Picard CollectHsMetrics FOLD_80_BASE_PENALTY."
        format: "{:.3f}"
        placement: 155
      Mean_Target_Coverage:
        title: "Mean Target Coverage"
        description: "Mean coverage over target intervals (x). Source: Picard CollectHsMetrics MEAN_TARGET_COVERAGE."
        format: "{:.2f}"
        placement: 165
      PCT_Selected_Bases:
        title: "Proportion Selected Bases"
        description: "Proportion of aligned bases mapping on or near target intervals (0-1). Source: Picard CollectHsMetrics PCT_SELECTED_BASES."
        format: "{:.4f}"
        placement: 170
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

    // Set defaults for optional params
    if (!params.architecture) params.architecture = "x86"

    // Coverage-targeted subsampling (Ultima/CRAM only). When set, every CRAM
    // biosample is subsampled to this single target coverage and renamed
    // `{biosample}_cov{NN}x`. Requires a `custom_mean_coverage` column in
    // input_csv. One target per pipeline run keeps each launch's cost
    // structure identical to a standard basej-wgs run; MCP queues N
    // launches when N coverage tiers are needed.
    def coverage_target = null
    if (params.coverage_target != null && params.coverage_target.toString().trim() != "" &&
            params.coverage_target.toString().trim() != "0") {
        try {
            coverage_target = (params.coverage_target as Number).intValue()
        } catch (Exception e) {
            exit 1, "ERROR: params.coverage_target must be a positive integer; got '${params.coverage_target}'"
        }
        if (coverage_target <= 0) {
            exit 1, "ERROR: params.coverage_target must be a positive integer; got ${params.coverage_target}"
        }
    }
    def use_coverage_subsample = (coverage_target != null)
    if (use_coverage_subsample) {
        log.info "Coverage subsampling enabled: target=${coverage_target}x. Each CRAM biosample will be renamed to {biosample}_cov${String.format('%02d', coverage_target)}x and subsampled."
    }

    // Parse input CSV with smart branch detection: CRAM or FASTQ
    ch_raw_csv = channel.fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)

    // Branch: detect CRAM vs FASTQ based on column content
    ch_split = ch_raw_csv.branch {
        cram: it.cram && it.cram.trim() != ""
        fastq: it.read1 && it.read2 && it.read1.trim() != "" && it.read2.trim() != ""
    }

    // CRAM path: biosampleName + CRAM + auto-derived CRAI + read_type="CRAM"
    // When coverage subsampling is enabled, validate custom_mean_coverage
    // and emit one record per CRAM with renamed `{biosample}_cov{NN}x`
    // and a per-record proportion = coverage_target / custom_mean_coverage.
    // The proportion is split into a separate sample-keyed channel so the
    // standard 4-tuple shape is preserved downstream.
    if (use_coverage_subsample) {
        ch_cram_with_meta = ch_split.cram
            .map { row ->
                if (!row.custom_mean_coverage || !row.custom_mean_coverage.toString().trim() ||
                        !row.custom_mean_coverage.toString().trim().isFloat()) {
                    exit 1, "ERROR: coverage_target is set but biosample '${row.biosampleName}' " +
                            "is missing a numeric 'custom_mean_coverage' column in input_csv. " +
                            "Add a 'custom_mean_coverage' value (e.g. 30.5) for every CRAM row, or unset coverage_target."
                }
                def mean_cov  = row.custom_mean_coverage.toString().trim() as float
                def proportion
                def cov_label
                if (coverage_target >= mean_cov) {
                    // Cannot upsample: this biosample's native coverage is at or
                    // below the requested target. Rather than aborting the whole
                    // run, pass it through at native coverage (proportion = 1.0,
                    // no subsampling) and label it with its actual (floored)
                    // coverage so downstream analysis does not mistake it for a
                    // true ${coverage_target}x sample.
                    proportion = 1.0f
                    cov_label  = mean_cov.intValue()
                    log.warn "Biosample '${row.biosampleName}': custom_mean_coverage (${mean_cov}x) <= coverage_target (${coverage_target}x); cannot upsample. " +
                             "Passing through at native coverage and labeling _cov${String.format('%02d', cov_label)}x (not _cov${String.format('%02d', coverage_target)}x)."
                } else {
                    proportion = (coverage_target as float) / mean_cov
                    cov_label  = coverage_target
                }
                def fanout_name = "${row.biosampleName}_cov${String.format('%02d', cov_label)}x"
                def cram_file   = file(row.cram)
                def crai_file   = file(row.cram + '.crai')
                // 5-tuple: carries proportion; split right after.
                [fanout_name, cram_file, crai_file, "CRAM", proportion as float]
            }

        // Standard 4-tuple shape for consistency with FASTQ branch
        ch_cram_reads = ch_cram_with_meta
            .map { fanout_name, cram, crai, read_type, prop -> [fanout_name, cram, crai, read_type] }
        // Keyed proportions for SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION
        ch_cram_proportions = ch_cram_with_meta
            .map { fanout_name, cram, crai, read_type, prop -> [fanout_name, prop] }
    } else {
        ch_cram_reads = ch_split.cram
            .map { row ->
                def cram_file = file(row.cram)
                def crai_file = file(row.cram + '.crai')
                [row.biosampleName, cram_file, crai_file, "CRAM"]
            }
        ch_cram_proportions = channel.empty()
    }

    // FASTQ path: multi-lane support via pipe-delimited ("|") read1/read2 paths
    //   Single-lane: biosampleName,s3://bucket/R1.fastq.gz,s3://bucket/R2.fastq.gz
    //   Multi-lane:  biosampleName,s3://bucket/L001_R1.fq.gz|s3://bucket/L002_R1.fq.gz,s3://bucket/L001_R2.fq.gz|s3://bucket/L002_R2.fq.gz
    ch_fastq_branched = ch_split.fastq
        .branch {
            multilane: it.read1.contains('|')
            singlelane: true
        }

    // Single-lane: one R1, one R2 — pass through directly with read_type
    ch_fastq_singlelane = ch_fastq_branched.singlelane
        .map { row -> [row.biosampleName, [file(row.read1), file(row.read2)], "FASTQ"] }

    // Multi-lane: split on "|", collect file objects, flatten for cat
    MERGE_MULTILANE_FASTQ(
        ch_fastq_branched.multilane
            .map { row ->
                def r1_files = row.read1.tokenize('|').collect { file(it.trim()) }
                def r2_files = row.read2.tokenize('|').collect { file(it.trim()) }
                [row.biosampleName, r1_files + r2_files]
            }
    )

    ch_fastq_reads = ch_fastq_singlelane
        .mix(MERGE_MULTILANE_FASTQ.out.reads.map { sample_id, reads -> [sample_id, reads, "FASTQ"] })

    // Merged channel: [sample_name, reads, read_type]
    // Note: For CRAM, reads = [cram_file]; for FASTQ, reads = [r1, r2]
    ch_reads_merged = ch_cram_reads.mix(ch_fastq_reads)

    ch_reads_merged.view { sample -> "Processing sample: ${sample[0]}, read_type: ${sample[sample.size()-1]}" }
    ch_reads_merged.ifEmpty { exit 1, "ERROR: No valid CRAM or FASTQ samples found in --input_csv" }

    // Mode prefix for output naming: wgs → wgsqc_*, exome → wesqc_*
    def mode_prefix = params.mode == 'exome' ? 'wes' : 'wgs'

    // Resolve mode-specific target intervals
    def wgs_or_target_intervals
    if (params.mode == 'exome') {
        wgs_or_target_intervals = params.genomes[params.genome][params.exome_panel]['wgs_or_target_intervals']
    } else {
        wgs_or_target_intervals = params.genomes[params.genome]['wgs_or_target_intervals']
    }

    // Step 1: Optional subsampling
    //   FASTQ:   SEQKIT_SAMPLE (read-count cap, max_total_reads)
    //   CRAM:    SAMTOOLS_SUBSAMPLE_CRAM (read-count cap)
    //              OR
    //            SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION (when coverage_target is
    //            set — proportion derived from custom_mean_coverage in CSV)
    if (use_coverage_subsample) {
        // Coverage-targeted subsampling: ignore skip_subsampling/max_total_reads
        // for the CRAM branch (FASTQ branch is untouched). Route every CRAM
        // record through the proportion-based subsampler.
        ch_split_for_subsample = ch_reads_merged.branch {
            cram: it[it.size()-1] == "CRAM"
            fastq: it[it.size()-1] == "FASTQ"
        }

        // Join records with their proportions (keyed by fanout_name)
        ch_cram_for_proportion = ch_split_for_subsample.cram
            .map { sample, cram, crai, read_type -> [sample, cram, crai] }
            .join(ch_cram_proportions)
            // [sample, cram, crai, proportion]

        SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION(
            ch_cram_for_proportion,
            params.samtools_seed
        )

        // FASTQ branch follows the existing skip_subsampling logic
        if (!params.skip_subsampling) {
            SEQKIT_SAMPLE(
                ch_split_for_subsample.fastq.map { sample, reads, read_type -> [sample, reads, params.max_total_reads] },
                params.seqkit_sample_seed
            )
            ch_aligned_reads = SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION.out.reads
                .map { sample, cram, crai -> [sample, cram, "CRAM"] }
                .mix(SEQKIT_SAMPLE.out.reads.map { sample, reads -> [sample, reads, "FASTQ"] })
        } else {
            ch_aligned_reads = SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION.out.reads
                .map { sample, cram, crai -> [sample, cram, "CRAM"] }
                .mix(ch_split_for_subsample.fastq.map { items -> [items[0], items[1], items[2]] })
        }

        // Coverage subsampling doesn't write read_count files; placeholder
        // so the parquet step's path() input is non-empty.
        ch_readcount_metrics = channel.of(file("${projectDir}/assets/dummy_file.txt")).collect()
    } else if (!params.skip_subsampling) {
        ch_split_for_subsample = ch_reads_merged.branch {
            cram: it[it.size()-1] == "CRAM"  // read_type is last element
            fastq: it[it.size()-1] == "FASTQ"
        }

        // CRAM subsampling: extra inputs required (cram, crai files explicitly)
        SAMTOOLS_SUBSAMPLE_CRAM(
            ch_split_for_subsample.cram.map { sample, cram, crai, read_type ->
                [sample, cram, crai, params.max_total_reads]
            },
            params.samtools_seed
        )

        // FASTQ subsampling: use standard FASTQ format [sample, [r1, r2]]
        SEQKIT_SAMPLE(
            ch_split_for_subsample.fastq.map { sample, reads, read_type -> [sample, reads, params.max_total_reads] },
            params.seqkit_sample_seed
        )

        // Merge subsampled outputs back to common format [sample, reads, read_type]
        ch_aligned_reads = SAMTOOLS_SUBSAMPLE_CRAM.out.reads
            .map { sample, cram, crai -> [sample, cram, "CRAM"] }
            .mix(SEQKIT_SAMPLE.out.reads.map { sample, reads -> [sample, reads, "FASTQ"] })
        
        // Collect read count files for Ultima/CRAM only
        // Use empty dummy file as fallback if no CRAM files (e.g., Illumina FASTQ-only run)
        ch_readcount_metrics = SAMTOOLS_SUBSAMPLE_CRAM.out.read_counts_file
            .collect()
            .ifEmpty([file("${projectDir}/assets/dummy_file.txt")])
    } else {
        // No subsampling: use reads as-is but normalize tuple arity to [sample, reads, read_type]
        // CRAM records come as [sample, cram, crai, "CRAM"], drop crai to match SENTIEON_ALIGN_DEDUP input spec
        ch_aligned_reads = ch_reads_merged.map { items ->
            if (items[items.size()-1] == "CRAM") {
                // CRAM: [sample, cram, crai, "CRAM"] → [sample, cram, "CRAM"]
                [items[0], items[1], items[3]]
            } else {
                // FASTQ: [sample, [r1, r2], "FASTQ"] → [sample, [r1, r2], "FASTQ"] (no change)
                items
            }
        }
        // Use empty dummy file as a placeholder to avoid empty channel error in Nextflow path() input
        // Python glob for *_read_counts.txt won't match it, so it's effectively ignored
        ch_readcount_metrics = channel.of(file("${projectDir}/assets/dummy_file.txt")).collect()
    }

    // Step 2: Branch by read_type
    //   FASTQ → SENTIEON_ALIGN_DEDUP → BAM → SENTIEON_DRIVER_METRICS
    //   CRAM  → SENTIEON_METRICS_CRAM (skip dedup; Ultima `demux` pre-marks dups)
    ch_aligned_branched = ch_aligned_reads.branch {
        cram:  it[2] == "CRAM"
        fastq: it[2] == "FASTQ"
    }

    // CRAM path needs the .crai too — re-attach from upstream channels
    if (use_coverage_subsample) {
        ch_cram_for_metrics = SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION.out.reads  // [sample, cram, crai]
    } else if (!params.skip_subsampling) {
        ch_cram_for_metrics = SAMTOOLS_SUBSAMPLE_CRAM.out.reads  // [sample, cram, crai]
    } else {
        ch_cram_for_metrics = ch_reads_merged
            .filter { it[it.size()-1] == "CRAM" }
            .map { sample, cram, crai, read_type -> [sample, cram, crai] }
    }

    // Step 2/3: Alignment, dedup and QC metrics.
    //   Toggle via params.pipeline_tool: 'sentieon' (proprietary — SENTIEON_ALIGN_DEDUP,
    //   SENTIEON_METRICS_CRAM and SENTIEON_DRIVER_METRICS) or the open-source path
    //   (BWAMEM2_ALIGN_DEDUP_METRICS + PICARD_METRICS_CRAM). Open-source is x86-only.
    //   Neutral channel handles below feed the existing downstream wiring unchanged:
    //     ch_bam               — FASTQ-derived dedup BAM  [sample, bam, bai]
    //     ch_dedup_metrics     — dedup metrics tuple       [sample, dedup_txt]
    //     ch_dedup_metrics_flat / ch_driver_metrics_flat / ch_cram_metrics_flat
    //     ch_driver_metrics_tuple / ch_cram_metrics_tuple  — for publish:
    if (params.pipeline_tool == 'sentieon') {
        SENTIEON_METRICS_CRAM(
            ch_cram_for_metrics,
            params.reference,
            params.base_metrics_intervals,
            wgs_or_target_intervals,
            params.mode
        )

        // FASTQ path: alignment + dedup
        SENTIEON_ALIGN_DEDUP(
            ch_aligned_branched.fastq,
            params.reference,
            params.platform
        )

        // QC Metrics for FASTQ-derived BAMs
        SENTIEON_DRIVER_METRICS(
            SENTIEON_ALIGN_DEDUP.out.bam,
            params.reference,
            params.base_metrics_intervals,
            wgs_or_target_intervals,
            params.mode
        )

        ch_bam                  = SENTIEON_ALIGN_DEDUP.out.bam
        ch_dedup_metrics        = SENTIEON_ALIGN_DEDUP.out.dedup_metrics
        ch_dedup_metrics_flat   = SENTIEON_ALIGN_DEDUP.out.dedup_metrics_flat
        ch_driver_metrics_flat  = SENTIEON_DRIVER_METRICS.out.metrics_flat
        ch_cram_metrics_flat    = SENTIEON_METRICS_CRAM.out.metrics_flat
        ch_driver_metrics_tuple = SENTIEON_DRIVER_METRICS.out.metrics_tuple
        ch_cram_metrics_tuple   = SENTIEON_METRICS_CRAM.out.metrics_tuple
    } else {
        // Open-source path: BWA-MEM2 → samtools markdup → Picard suite, all in
        // one process for FASTQ. Uses the pre-built BWA-MEM2 index staged at
        // params.bwamem2_reference (built by the nf-bwamem2-index-module workflow).
        BWAMEM2_ALIGN_DEDUP_METRICS(
            ch_aligned_branched.fastq,
            params.bwamem2_reference,
            params.reference,
            params.base_metrics_intervals,
            wgs_or_target_intervals,
            params.mode,
            params.platform
        )

        // CRAM path: metrics-only on the pre-aligned CRAM (Picard reads CRAM via ref)
        PICARD_METRICS_CRAM(
            ch_cram_for_metrics,
            params.reference,
            params.base_metrics_intervals,
            wgs_or_target_intervals,
            params.mode
        )

        ch_bam                  = BWAMEM2_ALIGN_DEDUP_METRICS.out.bam
        ch_dedup_metrics        = BWAMEM2_ALIGN_DEDUP_METRICS.out.dedup_metrics
        // dedup_sentieonmetrics.txt is already captured by metrics_flat's
        // *sentieonmetrics* glob, so keep this empty to avoid a duplicate
        // input-file-name collision in WGS_QC_METRICS_TO_PARQUET.
        ch_dedup_metrics_flat   = Channel.empty()
        ch_driver_metrics_flat  = BWAMEM2_ALIGN_DEDUP_METRICS.out.metrics_flat
        ch_cram_metrics_flat    = PICARD_METRICS_CRAM.out.metrics_flat
        ch_driver_metrics_tuple = BWAMEM2_ALIGN_DEDUP_METRICS.out.metrics_tuple
        ch_cram_metrics_tuple   = PICARD_METRICS_CRAM.out.metrics_tuple
    }

    // Step 4: Picard HS metrics (exome only, FASTQ path only — Ultima exome
    // is rare and PE-oriented Picard HsMetrics on single-end CRAM is noisy)
    // Collect all metrics into a single channel for WGS_QC_METRICS_TO_PARQUET.
    // FASTQ-path dedup_metrics + driver metrics; CRAM-path metrics include
    // empty placeholder dedup_sentieonmetrics so the parser is uniform.
    ch_all_metrics = ch_dedup_metrics_flat
        .mix(ch_driver_metrics_flat)
        .mix(ch_cram_metrics_flat)

    if (params.mode == 'exome') {
        // Run Picard HsMetrics for both FASTQ-derived BAMs and Ultima CRAMs.
        // Picard CollectHsMetrics reads CRAM transparently when given a
        // reference sequence, so the same process handles both inputs.
        ch_for_hsmetrics = ch_bam
            .mix(ch_cram_for_metrics)
        PICARD_COLLECTHSMETRICS(
            ch_for_hsmetrics,
            params.reference,
            wgs_or_target_intervals
        )
        ch_all_metrics = ch_all_metrics.mix(PICARD_COLLECTHSMETRICS.out.metrics_flat)
    }

    // Step 5: Parse all metrics → Parquet + MultiQC TSVs
    // Pass read count files if available (Ultima/CRAM only), or empty channel
    WGS_QC_METRICS_TO_PARQUET(
        ch_all_metrics.collect(),
        ch_readcount_metrics,
        params.mode,
        mode_prefix,
        params.genome,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        params.pipeline_user
    )

    // Step 6: QC Plots + final parquets with QC status
    WGS_QC_PLOTS(
        WGS_QC_METRICS_TO_PARQUET.out.summary_tsv,
        params.mode,
        mode_prefix,
        params.workspace,
        params.workflow_id
    )

    // Step 7: MultiQC Report
    // Collect all file types together: raw metrics (for Picard/Sentieon parsing) +
    // processed outputs (custom QC tables/plots/scores from WGS_QC_PLOTS).
    // Raw metrics files stage with unique per-sample names, avoiding collisions.
    ch_mqc_inputs = ch_all_metrics
        .mix(WGS_QC_PLOTS.out.mqc_metrics)
        .mix(WGS_QC_PLOTS.out.plots_jpg)
        .mix(WGS_QC_PLOTS.out.scores_summary)
        .collect()

    MULTIQC(
        ch_mqc_inputs,
        params.mode,
        params.dataset_id,
        params.workspace,
        params.workflow_id,
        workflow.manifest.version,
        file("${projectDir}/assets/bioskryb_logo-tagline.png")
    )

    publish:
    bam_files = ch_bam
        .map { sample_name, bam, bai -> [biosampleName: sample_name, bam: bam, bai: bai] }

    // Subsampled CRAMs (coverage_target only) — published with renamed
    // biosample names like `{parent}_cov{NN}x` so MCP can queue downstream
    // DV/vcfeval runs per coverage tier. Empty channel when not subsampling
    // by coverage.
    cram_files = use_coverage_subsample
        ? SAMTOOLS_SUBSAMPLE_CRAM_PROPORTION.out.reads
            .map { sample_name, cram, crai -> [biosampleName: sample_name, cram: cram, crai: crai] }
        : channel.empty()

    dedup_metrics = ch_dedup_metrics
        .map { sample_name, metrics -> [biosampleName: sample_name, metrics: metrics] }

    sentieon_metrics = ch_driver_metrics_tuple
        .mix(ch_cram_metrics_tuple)
        .map { sample_name, metrics -> [biosampleName: sample_name, metrics: metrics] }

    wgsqc_summary     = WGS_QC_PLOTS.out.parquet
    wgsqc_summary_tsv = WGS_QC_PLOTS.out.summary_tsv
    wgsqc_plots       = channel.empty()
        .mix(WGS_QC_PLOTS.out.plots_pdf)
        .mix(WGS_QC_PLOTS.out.plots_jpg)
    wgsqc_scores      = channel.empty()
        .mix(WGS_QC_PLOTS.out.scores)
        .mix(WGS_QC_PLOTS.out.scores_summary)
    per_biosample_status = WGS_QC_PLOTS.out.per_biosample_status
    multiqc_report    = MULTIQC.out.report
}

// ============================================================================
// OUTPUT CONFIGURATION
// ============================================================================
output {
    bam_files {
        path "bam/${params.workspace}/dna/tool=${params.pipeline_tool == 'sentieon' ? 'sentieon-202503-02' : 'bwa-mem2'}/pipeline=wgsqc"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/bam.csv"
            header true
        }
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "bam",
             tool:         "${params.pipeline_tool == 'sentieon' ? 'sentieon' : 'bwa-mem2'}".toString(),
             reference:    params.genome
    }

    cram_files {
        path "cram/${params.workspace}/dna/tool=samtools-subsampled"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/cram.csv"
            header true
        }
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "cram",
             tool:         "samtools",
             reference:    params.genome
    }

    dedup_metrics {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/dedup"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/dedup_metrics.csv"
            header true
        }
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "dedup_metrics",
             tool:         "${params.pipeline_tool == 'sentieon' ? 'sentieon' : 'samtools'}".toString(),
             reference:    params.genome
    }

    sentieon_metrics {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/sentieon"
        index {
            path "workflow_outputs/${params.workspace}/${params.workflow_id}/index/sentieon_metrics.csv"
            header true
        }
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "qc_metrics",
             tool:         "${params.pipeline_tool == 'sentieon' ? 'sentieon' : 'picard'}".toString(),
             reference:    params.genome
    }

    wgsqc_summary {
        path "tables"
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_summary".toString(),
             reference:    params.genome
    }

    wgsqc_summary_tsv {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/${params.mode == 'exome' ? 'wes' : 'wgs'}qc_metrics"
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_summary_tsv".toString()
    }

    wgsqc_plots {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/qc_plots"
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_plots".toString()
    }

    wgsqc_scores {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/metrics/${params.mode == 'exome' ? 'wes' : 'wgs'}qc_metrics"
        tags workspace:    params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             molecule_type: "dna",
             artifact:     "${params.mode == 'exome' ? 'wes' : 'wgs'}qc_scores".toString()
    }

    multiqc_report {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/reports"
        tags workspace:   params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id: params.workflow_id,
             pipeline:    "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             artifact:    "multiqc_report"
    }

    per_biosample_status {
        path "workflow_outputs/${params.workspace}/${params.workflow_id}/index"
        tags workspace:   params.workspace,
             dataset_id:  params.dataset_id,
             workflow_id:  params.workflow_id,
             pipeline:     "${params.mode == 'exome' ? 'basej-wesqc' : 'basej-wgsqc'}".toString(),
             artifact:    "per_biosample_status"
    }
}
