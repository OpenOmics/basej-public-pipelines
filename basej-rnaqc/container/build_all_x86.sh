#!/usr/bin/env bash
# Build all x86 basej-rnaqc containers locally, one log per image + a summary.
set -u

BASE="$(cd "$(dirname "$0")" && pwd)"
LOGDIR="$BASE/_build_logs"
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/SUMMARY.log"
: > "$SUMMARY"

# folder:image_tag  (x86 only)
BUILDS=(
  "ubuntu:basejumper_ubuntu_24.04_stable"
  "seqkit:basejumper_seqkit-2.13.0"
  "custom_r_nf_rnaseq:basejumper_custom_r_nf_rnaseq_0.10"
  "custom_r_qcplots:basejumper_custom_r_qcplots_0.2.3"
  "gene_body_coverage:basejumper_gene_body_coverage_0.2.1"
)

echo "BUILD STARTED: $(date)" | tee -a "$SUMMARY"
for entry in "${BUILDS[@]}"; do
  folder="${entry%%:*}"
  tag="${entry##*:}"
  log="$LOGDIR/${folder}.log"
  echo "[$(date +%H:%M:%S)] BUILDING $folder -> $tag (log: $log)" | tee -a "$SUMMARY"
  if docker build -t "$tag" "$BASE/$folder" > "$log" 2>&1; then
    echo "[$(date +%H:%M:%S)] OK      $folder -> $tag" | tee -a "$SUMMARY"
  else
    echo "[$(date +%H:%M:%S)] FAILED  $folder -> $tag (see $log)" | tee -a "$SUMMARY"
  fi
done
echo "BUILD FINISHED: $(date)" | tee -a "$SUMMARY"
