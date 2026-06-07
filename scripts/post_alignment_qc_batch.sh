#!/bin/bash
# =============================================================================
# fastq2tracks v3.1.1 — Post-alignment QC module (deepTools-based, parallelised)
# Replaces ChIPQC with a robust, crash-proof deepTools + samtools/bedtools QC.
# Compatible assays: ChIP-seq, CUT&RUN, CUT&Tag, ChIPmentation, ATAC-seq
#
# Parallelism strategy (safe for shared servers):
# Phase 1: FRiP per sample → background jobs, throttled to QC_PARALLEL_JOBS
# Phase 2: bamCoverage (karyogram) → background jobs, throttled to QC_PARALLEL_JOBS
#          each bamCoverage uses QC_THREADS_PER_JOB threads (default 2)
# Phase 3: multiBamSummary bins runs alongside plotFingerprint (already was)
#          plotCorrelation pearson/spearman run in parallel after bins NPZ ready
#          multiBamSummary BED-file starts immediately after consensus built
# Phase 4: FRiP consensus → background jobs, throttled to QC_PARALLEL_JOBS
#
# New config.conf variables (all optional, have safe defaults):
# QC_PARALLEL_JOBS    — max background jobs at once (default: 4)
# QC_THREADS_PER_JOB  — threads per parallel bamCoverage job (default: 2)
# THREADS_DEEPTOOLS   — threads for serial deepTools steps (unchanged)
#
# Tuning guide (500 GB RAM / 140 logical cores server):
# Conservative (other users present):
#   QC_PARALLEL_JOBS=4, QC_THREADS_PER_JOB=2, THREADS_DEEPTOOLS=16
# Moderate (dedicated run):
#   QC_PARALLEL_JOBS=8, QC_THREADS_PER_JOB=4, THREADS_DEEPTOOLS=32
# Aggressive (server to yourself):
#   QC_PARALLEL_JOBS=12, QC_THREADS_PER_JOB=6, THREADS_DEEPTOOLS=48
#
# Usage (called by fastq2tracks.sh as Step 10):
#   bash scripts/post_alignment_qc_batch.sh \
#        <samplesheet.csv> <bam_dir> <peaks_dir> <bigwig_dir> <out_dir>
# =============================================================================
set -uo pipefail   # NOTE: -e intentionally omitted — errors handled per-step

# ── Config loading ─────────────────────────────────────────────────────────────
_load_config() {
  if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
    source "${F2T_CONFIG}"
  else
    local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
    local _c="${_d}/../config/config.conf"
    [[ -f "$_c" ]] && source "$_c" || {
      echo "ERROR: config.conf not found. Export F2T_CONFIG or pass --config." >&2
      exit 1
    }
  fi
}
_load_config

# ── Arguments ──────────────────────────────────────────────────────────────────
SAMPLESHEET="${1:?Usage: post_alignment_qc_batch.sh <samplesheet> <bam_dir> <peaks_dir> <bigwig_dir> <out_dir>}"
BAM_DIR="${2:?BAM_DIR required}"
PEAKS_DIR="${3:?PEAKS_DIR required}"
BIGWIG_DIR="${4:?BIGWIG_DIR required}"
OUT_DIR="${5:?OUT_DIR required}"

# ── Parallelism settings ───────────────────────────────────────────────────────
THREADS="${THREADS_DEEPTOOLS:-8}"
QC_PARALLEL_JOBS="${QC_PARALLEL_JOBS:-4}"
QC_THREADS_PER_JOB="${QC_THREADS_PER_JOB:-2}"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Output directory structure ─────────────────────────────────────────────────
QC_TABLES="${OUT_DIR}/tables"
QC_PLOTS="${OUT_DIR}/plots"
QC_CHRPLOTS="${OUT_DIR}/plots/chromosome_coverage"
QC_MATRICES="${OUT_DIR}/matrices"
QC_LOGS="${OUT_DIR}/logs"
QC_PEAKS="${OUT_DIR}/peak_sets"
QC_DT="${OUT_DIR}/deeptools"
QC_TMP="${OUT_DIR}/.tmp_frip"
mkdir -p "$QC_TABLES" "$QC_PLOTS" "$QC_CHRPLOTS" "$QC_MATRICES" \
         "$QC_LOGS" "$QC_PEAKS" "$QC_DT" "$QC_TMP"

MAIN_LOG="${QC_LOGS}/post_qc_${TIMESTAMP}.log"
SUMMARY_TSV="${QC_TABLES}/qc_summary.tsv"
WARNINGS_TSV="${QC_TABLES}/qc_warnings.tsv"

# ── Logging ────────────────────────────────────────────────────────────────────
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $1" | tee -a "$MAIN_LOG"; }

log "=== Post-alignment QC module started (v3.1.1 — parallelised) ==="
log "Samplesheet      : $SAMPLESHEET"
log "BAM dir          : $BAM_DIR"
log "Peaks dir        : $PEAKS_DIR"
log "BigWig dir       : $BIGWIG_DIR"
log "Output dir       : $OUT_DIR"
log "THREADS_DEEPTOOLS: $THREADS (serial deepTools steps)"
log "QC_PARALLEL_JOBS : $QC_PARALLEL_JOBS (background job slots)"
log "QC_THREADS_PER_JOB: $QC_THREADS_PER_JOB (threads per bamCoverage job)"
log "Peak bamCoverage CPU ceiling: $((QC_PARALLEL_JOBS * QC_THREADS_PER_JOB)) threads"

# ── Tool checks ────────────────────────────────────────────────────────────────
MISSING_TOOLS=()
for tool in samtools bedtools bamCoverage multiBamSummary multiBigwigSummary \
            plotCorrelation plotPCA plotFingerprint computeMatrix \
            plotHeatmap plotProfile; do
  command -v "$tool" &>/dev/null || MISSING_TOOLS+=("$tool")
done
[[ ${#MISSING_TOOLS[@]} -gt 0 ]] && \
  warn "Missing tools (some steps may skip): ${MISSING_TOOLS[*]}"

KARYOGRAM_PY="${SCRIPT_DIR}/plot_chrom_coverage.py"
[[ ! -f "$KARYOGRAM_PY" ]] && warn "plot_chrom_coverage.py not found — karyogram plots will be skipped" && KARYOGRAM_PY=""
if [[ -n "$KARYOGRAM_PY" ]]; then
  python3 -c "import matplotlib, numpy" 2>/dev/null \
    || { warn "matplotlib/numpy not available — karyogram plots will be skipped"; KARYOGRAM_PY=""; }
fi

# ── TSV headers ────────────────────────────────────────────────────────────────
echo -e "sample_id\tassay\tbam_file\ttotal_reads\tduplication_pct\tbl_filtered_reads\tmito_reads\tmito_pct\tn_narrow_peaks\tn_broad_peaks\tpeak_width_median\tpeak_max_signal\tfrip_narrow\tfrip_broad\tqc_status\tnotes" \
  > "$SUMMARY_TSV"
echo -e "sample_id\twarning_type\tvalue\tthreshold\tmessage" > "$WARNINGS_TSV"

# ── Helper: thread-safe warning ────────────────────────────────────────────────
add_warning() {
  local sid="$1" wtype="$2" val="$3" thr="$4" msg="$5"
  echo -e "${sid}\t${wtype}\t${val}\t${thr}\t${msg}" >> "$WARNINGS_TSV"
  warn "$sid — $wtype: $msg"
}

# ── Helper: mitochondrial fraction ────────────────────────────────────────────
compute_mito_fraction() {
  local key="$1" fallback_bam="$2"
  local bam_mito="${BAM_DIR}/../dedupBams/${key}_dedup.bam"
  [[ ! -f "$bam_mito" ]] && bam_mito="$fallback_bam"
  if [[ ! -f "$bam_mito" ]]; then echo "NA NA"; return; fi
  local mito total
  mito=$(samtools idxstats "$bam_mito" 2>/dev/null \
    | awk '$1=="chrM"||$1=="MT"||$1=="chrMT"{sum+=$3}END{print sum+0}')
  total=$(samtools idxstats "$bam_mito" 2>/dev/null \
    | awk '{sum+=$3}END{print sum+0}')
  [[ "$total" -eq 0 ]] && echo "NA NA" && return
  awk "BEGIN{printf \"%d %.2f\", $mito, ($mito/$total)*100}"
}

# ── Helper: throttled background job runner ────────────────────────────────────
declare -a _PIDS=()
_reap_finished() {
  local new_pids=()
  for p in "${_PIDS[@]}"; do
    kill -0 "$p" 2>/dev/null && new_pids+=("$p")
  done
  _PIDS=("${new_pids[@]+${new_pids[@]}}")
}
_wait_slot() {
  local max="$1"
  while true; do
    _reap_finished
    [[ "${#_PIDS[@]}" -lt "$max" ]] && break
    sleep 1
  done
}
_drain() {
  for p in "${_PIDS[@]+${_PIDS[@]}}"; do
    wait "$p" 2>/dev/null || true
  done
  _PIDS=()
}

# =============================================================================
# Phase 1: Per-sample metrics (reads, dups, mito, peaks)
# FRiP computed in parallel background jobs — throttled to QC_PARALLEL_JOBS
# =============================================================================
log "=== Phase 1: Per-sample metrics ==="

declare -A seen_keys
declare -a SAMPLE_KEYS SAMPLE_BAMS SAMPLE_NARROW SAMPLE_BROAD SAMPLE_BIGWIGS
declare -A SAMPLE_GENOME SAMPLE_ASSAY
declare -A SAMPLE_READS SAMPLE_DUP SAMPLE_MITO_R SAMPLE_MITO_P
declare -A SAMPLE_N_NARROW SAMPLE_N_BROAD SAMPLE_PEAK_WIDTH SAMPLE_PEAK_SIG
declare -A SAMPLE_QC_STATUS

while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment \
      cell_type rep tech_rep is_ctrl ctrl_id macs2_mode blacklist rest; do

  [[ "$sid" == "sample_id" ]] && continue
  sid="${sid//\"/}"; rep="${rep//\"/}"; is_ctrl="${is_ctrl//\"/}"
  assay="${assay//\"/}"; genome="${genome//\"/}"
  [[ "${is_ctrl,,}" == "true" || "$is_ctrl" == "1" ]] && continue

  KEY="${sid}_bioR${rep}"
  [[ -n "${seen_keys[$KEY]+x}" ]] && continue
  seen_keys["$KEY"]=1

  BAM="${BAM_DIR}/${KEY}_dedup_blFilt.bam"
  [[ ! -f "$BAM" ]] && BAM="${BAM_DIR}/${KEY}_dedup.bam"
  if [[ ! -f "$BAM" ]]; then
    warn "$KEY: BAM not found — recording as FAILED"
    echo -e "${KEY}\t${assay}\tNA\t0\tNA\t0\tNA\tNA\t0\t0\tNA\tNA\tNA\tNA\tFAILED\tBAM not found" >> "$SUMMARY_TSV"
    continue
  fi

  log "  Collecting stats: $KEY (assay=$assay genome=$genome)"
  SAMPLE_GENOME["$KEY"]="$genome"
  SAMPLE_ASSAY["$KEY"]="${assay,,}"

  DUP_METRICS=$(find \
    "${OUT_DIR}/../logs/picard" \
    "${BAM_DIR}/../logs/picard" \
    "${BAM_DIR}/../dedupBams" \
    -name "${KEY}_dup_metrics.txt" 2>/dev/null | head -1)
  DUP_PCT="NA"
  [[ -f "$DUP_METRICS" ]] && \
    DUP_PCT=$(awk '/^LIBRARY/{found=1;next} found&&NF>0{printf "%.2f",$9*100;exit}' \
      "$DUP_METRICS" 2>/dev/null || echo "NA")

  TOTAL_READS=$(samtools view -c -F 4 "$BAM" 2>/dev/null || echo 0)
  SAMPLE_READS["$KEY"]="$TOTAL_READS"
  SAMPLE_DUP["$KEY"]="$DUP_PCT"

  [[ "$TOTAL_READS" -lt 1000000 ]] && \
    add_warning "$KEY" "LOW_READS" "$TOTAL_READS" "1000000" "Very low aligned reads (<1M)"
  if [[ "$DUP_PCT" != "NA" ]]; then
    dup_int=${DUP_PCT%.*}
    [[ "$dup_int" -gt 80 ]] && \
      add_warning "$KEY" "HIGH_DUPLICATION" "$DUP_PCT" "80" "Duplication >80%"
  fi

  MITO_INFO=($(compute_mito_fraction "$KEY" "$BAM"))
  SAMPLE_MITO_R["$KEY"]="${MITO_INFO[0]:-NA}"
  SAMPLE_MITO_P["$KEY"]="${MITO_INFO[1]:-NA}"
  MITO_PCT="${MITO_INFO[1]:-NA}"
  if [[ "${assay,,}" =~ ^atac && "$MITO_PCT" != "NA" ]]; then
    mito_int=${MITO_PCT%.*}
    [[ "$mito_int" -ge 25 ]] && \
      add_warning "$KEY" "HIGH_MITO_ATAC" "$MITO_PCT" "25" \
        "ATAC-seq mito reads ≥25% — poor nuclear accessibility or degraded sample"
    [[ "$mito_int" -ge 10 && "$mito_int" -lt 25 ]] && \
      add_warning "$KEY" "ELEVATED_MITO_ATAC" "$MITO_PCT" "10" \
        "ATAC-seq mito reads ≥10% — check sample quality"
  fi

  NARROW_PEAK="${PEAKS_DIR}/per_replicate/${KEY}/narrow/${KEY}_peaks.narrowPeak"
  BROAD_PEAK="${PEAKS_DIR}/per_replicate/${KEY}/broad/${KEY}_peaks.broadPeak"
  N_NARROW=0; N_BROAD=0; PEAK_WIDTH_MED="NA"; PEAK_MAX_SIG="NA"

  if [[ -f "$NARROW_PEAK" && -s "$NARROW_PEAK" ]]; then
    N_NARROW=$(wc -l < "$NARROW_PEAK" 2>/dev/null || echo 0)
    PEAK_WIDTH_MED=$(awk '{print $3-$2}' "$NARROW_PEAK" | sort -n \
      | awk '{a[NR]=$0}END{print (NR%2==0)?(a[NR/2]+a[NR/2+1])/2:a[int(NR/2)+1]}' \
      2>/dev/null || echo "NA")
    PEAK_MAX_SIG=$(awk 'BEGIN{max=0}{if($7>max)max=$7}END{print max}' \
      "$NARROW_PEAK" 2>/dev/null || echo "NA")
  fi
  [[ -f "$BROAD_PEAK" && -s "$BROAD_PEAK" ]] && \
    N_BROAD=$(wc -l < "$BROAD_PEAK" 2>/dev/null || echo 0)

  SAMPLE_N_NARROW["$KEY"]="$N_NARROW"
  SAMPLE_N_BROAD["$KEY"]="$N_BROAD"
  SAMPLE_PEAK_WIDTH["$KEY"]="$PEAK_WIDTH_MED"
  SAMPLE_PEAK_SIG["$KEY"]="$PEAK_MAX_SIG"

  if [[ "$N_NARROW" -eq 0 && "$N_BROAD" -eq 0 ]]; then
    add_warning "$KEY" "ZERO_PEAKS" "0" "1" "No peaks called (narrow or broad)"
    SAMPLE_QC_STATUS["$KEY"]="NO_PEAKS"
  elif [[ "$N_NARROW" -lt 200 && "$N_BROAD" -lt 200 ]]; then
    add_warning "$KEY" "FEW_PEAKS" "$N_NARROW" "200" "Very few peaks (<200)"
    SAMPLE_QC_STATUS["$KEY"]="POOR"
  else
    SAMPLE_QC_STATUS["$KEY"]="PASS"
  fi

  BW="${BIGWIG_DIR}/${KEY}_dedup_blFilt_Snorm.bw"
  [[ ! -f "$BW" ]] && BW="${BIGWIG_DIR}/${KEY}_Snorm.bw"

  SAMPLE_KEYS+=("$KEY")
  SAMPLE_BAMS+=("$BAM")
  SAMPLE_NARROW+=("$NARROW_PEAK")
  SAMPLE_BROAD+=("$BROAD_PEAK")
  SAMPLE_BIGWIGS+=("$BW")

done < <(tail -n +2 "$SAMPLESHEET")

N_SAMPLES=${#SAMPLE_KEYS[@]}
log "Collected $N_SAMPLES IP samples"

if [[ $N_SAMPLES -eq 0 ]]; then
  warn "No valid samples — skipping deepTools cross-sample steps"
  exit 0
fi

# ── Pass 2: FRiP — parallel background jobs ────────────────────────────────────
log "  Computing FRiP in parallel (max $QC_PARALLEL_JOBS jobs)..."
_PIDS=()

for i in "${!SAMPLE_KEYS[@]}"; do
  KEY="${SAMPLE_KEYS[$i]}"
  BAM="${SAMPLE_BAMS[$i]}"
  NARROW_PEAK="${SAMPLE_NARROW[$i]}"
  BROAD_PEAK="${SAMPLE_BROAD[$i]}"
  TMP_OUT="${QC_TMP}/frip_${KEY}.tmp"
  [[ -f "$TMP_OUT" ]] && continue

  _wait_slot "$QC_PARALLEL_JOBS"

  (
    _frip() {
      local bam="$1" peak="$2"
      [[ ! -f "$bam" || ! -f "$peak" || ! -s "$peak" ]] && echo "NA" && return
      local total rip
      total=$(samtools view -c -F 4 "$bam" 2>/dev/null || echo 0)
      [[ "$total" -eq 0 ]] && echo "NA" && return
      rip=$(bedtools sort -i "$peak" 2>/dev/null \
        | bedtools merge -i stdin 2>/dev/null \
        | bedtools intersect -u -a "$bam" -b stdin -ubam 2>/dev/null \
        | samtools view -c 2>/dev/null || echo 0)
      awk "BEGIN {printf \"%.4f\", $rip / $total}"
    }
    fn=$(_frip "$BAM" "$NARROW_PEAK")
    fb=$(_frip "$BAM" "$BROAD_PEAK")
    printf '%s\t%s\t%s\n' "$KEY" "$fn" "$fb" > "$TMP_OUT"
  ) &
  _PIDS+=("$!")
done

log "  Waiting for FRiP jobs to finish..."
_drain
log "  FRiP jobs complete."

# ── Collect FRiP results + write summary TSV ───────────────────────────────────
declare -A FRIP_NARROW_MAP FRIP_BROAD_MAP
for f in "${QC_TMP}"/frip_*.tmp; do
  [[ -f "$f" ]] || continue
  read -r key fn fb < "$f"
  FRIP_NARROW_MAP["$key"]="${fn:-NA}"
  FRIP_BROAD_MAP["$key"]="${fb:-NA}"
done

for KEY in "${SAMPLE_KEYS[@]}"; do
  BAM_F="${SAMPLE_BAMS[0]}"
  for i in "${!SAMPLE_KEYS[@]}"; do
    [[ "${SAMPLE_KEYS[$i]}" == "$KEY" ]] && BAM_F="${SAMPLE_BAMS[$i]}" && break
  done
  FN="${FRIP_NARROW_MAP[$KEY]:-NA}"
  FB="${FRIP_BROAD_MAP[$KEY]:-NA}"

  if [[ "$FN" != "NA" ]]; then
    frip_pct=$(awk "BEGIN{printf \"%.0f\", $FN * 100}")
    [[ "$frip_pct" -lt 1 ]] && \
      add_warning "$KEY" "LOW_FRIP_NARROW" "$FN" "0.01" "FRiP (narrow) <1%"
  fi

  ASSAY="${SAMPLE_ASSAY[$KEY]:-unknown}"
  READS="${SAMPLE_READS[$KEY]:-0}"
  DUP="${SAMPLE_DUP[$KEY]:-NA}"
  MITO_R="${SAMPLE_MITO_R[$KEY]:-NA}"
  MITO_P="${SAMPLE_MITO_P[$KEY]:-NA}"
  N_NR="${SAMPLE_N_NARROW[$KEY]:-0}"
  N_BR="${SAMPLE_N_BROAD[$KEY]:-0}"
  PW="${SAMPLE_PEAK_WIDTH[$KEY]:-NA}"
  PS="${SAMPLE_PEAK_SIG[$KEY]:-NA}"
  QS="${SAMPLE_QC_STATUS[$KEY]:-UNKNOWN}"

  echo -e "${KEY}\t${ASSAY}\t${BAM_F}\t${READS}\t${DUP}\t${READS}\t${MITO_R}\t${MITO_P}\t${N_NR}\t${N_BR}\t${PW}\t${PS}\t${FN}\t${FB}\t${QS}\t" \
    >> "$SUMMARY_TSV"
  log "  $KEY: reads=$READS dup=$DUP mito=${MITO_P}% narrow=$N_NR broad=$N_BR frip_n=$FN"
done

rm -f "${QC_TMP}"/frip_*.tmp
log "Phase 1 complete."

# =============================================================================
# Phase 2a: Consensus peak set (fast bedtools merge)
# =============================================================================
log "=== Phase 2a: Consensus peak set (fast bedtools merge) ==="

NARROW_WITH_PEAKS=(); BROAD_WITH_PEAKS=()
for pk in "${SAMPLE_NARROW[@]}"; do [[ -f "$pk" && -s "$pk" ]] && NARROW_WITH_PEAKS+=("$pk"); done
for pk in "${SAMPLE_BROAD[@]}";  do [[ -f "$pk" && -s "$pk" ]] && BROAD_WITH_PEAKS+=("$pk");  done

MERGED_NARROW="${QC_PEAKS}/merged_narrow.bed"
MERGED_BROAD="${QC_PEAKS}/merged_broad.bed"
CONSENSUS_PEAK="${QC_PEAKS}/consensus_peaks.bed"

[[ ${#NARROW_WITH_PEAKS[@]} -gt 0 ]] && \
  cat "${NARROW_WITH_PEAKS[@]}" | cut -f1-3 | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin > "$MERGED_NARROW" 2>>"$MAIN_LOG" \
  && log "  Merged narrow: $(wc -l < "$MERGED_NARROW") regions"

[[ ${#BROAD_WITH_PEAKS[@]} -gt 0 ]] && \
  cat "${BROAD_WITH_PEAKS[@]}" | cut -f1-3 | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin > "$MERGED_BROAD" 2>>"$MAIN_LOG" \
  && log "  Merged broad: $(wc -l < "$MERGED_BROAD") regions"

PEAK_FILES_FOR_CONSENSUS=()
[[ -f "$MERGED_NARROW" && -s "$MERGED_NARROW" ]] && PEAK_FILES_FOR_CONSENSUS+=("$MERGED_NARROW")
[[ -f "$MERGED_BROAD"  && -s "$MERGED_BROAD"  ]] && PEAK_FILES_FOR_CONSENSUS+=("$MERGED_BROAD")

if [[ ${#PEAK_FILES_FOR_CONSENSUS[@]} -gt 0 ]]; then
  cat "${PEAK_FILES_FOR_CONSENSUS[@]}" | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin > "$CONSENSUS_PEAK" 2>>"$MAIN_LOG" \
  && log "  Consensus peaks: $(wc -l < "$CONSENSUS_PEAK") regions"
fi

# =============================================================================
# Phase 2b: Karyogram bedGraphs — parallel bamCoverage jobs
# =============================================================================
log "=== Phase 2b: Karyogram bedGraphs (parallel bamCoverage, max $QC_PARALLEL_JOBS jobs) ==="

KARYOGRAM_BG_LIST=()
_PIDS=()

for i in "${!SAMPLE_KEYS[@]}"; do
  KEY="${SAMPLE_KEYS[$i]}"
  BAM="${SAMPLE_BAMS[$i]}"
  GENOME="${SAMPLE_GENOME[$KEY]:-hg38}"
  [[ "${GENOME,,}" == "hg38" ]] && CHROM_SIZES="${CHROM_SIZES_HUMAN:-}" \
                                 || CHROM_SIZES="${CHROM_SIZES_MOUSE:-}"
  BG_OUT="${QC_CHRPLOTS}/${KEY}_100kb.bedGraph"
  CHR_TSV="${QC_CHRPLOTS}/${KEY}_per_chrom.tsv"

  if [[ ! -f "$CHR_TSV" ]]; then
    {
      echo -e "chrom\tchrom_reads\ttotal_reads\treads_per_mb"
      total=$(samtools view -c -F 4 "$BAM" 2>/dev/null || echo 1)
      samtools idxstats "$BAM" 2>/dev/null \
        | awk -v tot="$total" \
          '$1~/^chr[0-9XY]+$/ || $1=="chrM" || $1=="MT" {
             rpm=($2>0)?($3/$2*1e6):0
             printf "%s\t%d\t%d\t%.4f\n",$1,$3,tot,rpm
           }'
    } > "$CHR_TSV" 2>/dev/null
  fi

  if [[ -f "$BG_OUT" && -s "$BG_OUT" ]]; then
    log "  Skipping bamCoverage (exists): $KEY"
    KARYOGRAM_BG_LIST+=("$KEY:$BG_OUT")
    continue
  fi

  _wait_slot "$QC_PARALLEL_JOBS"

  (
    log "  bamCoverage 100kb start: $KEY"
    if bamCoverage \
        -b "$BAM" \
        --binSize 100000 \
        --normalizeUsing RPKM \
        --skipNonCoveredRegions \
        --outFileFormat bedgraph \
        -p "$QC_THREADS_PER_JOB" \
        -o "$BG_OUT" \
        >> "$MAIN_LOG" 2>&1; then
      log "  bamCoverage OK: $KEY"
    else
      warn "  bamCoverage FAILED: $KEY — trying bedtools fallback"
      if [[ -n "${CHROM_SIZES:-}" && -f "${CHROM_SIZES:-}" ]]; then
        local_tmp="${QC_CHRPLOTS}/tmp_tiles_${KEY}.bed"
        bedtools makewindows -g "$CHROM_SIZES" -w 100000 \
          | awk '$1~/^chr[0-9XY]+$/ || $1=="chrM"' \
          > "$local_tmp" 2>/dev/null
        bedtools coverage -a "$local_tmp" -b "$BAM" -counts \
          | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
          > "$BG_OUT" 2>/dev/null \
          && log "  bedtools fallback OK: $KEY" \
          || warn "  bedtools fallback also FAILED: $KEY"
        rm -f "$local_tmp"
      fi
    fi
  ) &
  _PIDS+=("$!")
done

log "  Waiting for all bamCoverage jobs..."
_drain
log "  bamCoverage jobs complete."

for i in "${!SAMPLE_KEYS[@]}"; do
  KEY="${SAMPLE_KEYS[$i]}"
  BG_OUT="${QC_CHRPLOTS}/${KEY}_100kb.bedGraph"
  [[ -f "$BG_OUT" && -s "$BG_OUT" ]] && KARYOGRAM_BG_LIST+=("$KEY:$BG_OUT")
done

CHR_SUMMARY_TSV="${QC_TABLES}/per_chromosome_reads.tsv"
{
  echo -e "sample_id\tchrom\tchrom_reads\ttotal_reads\treads_per_mb"
  for KEY in "${SAMPLE_KEYS[@]}"; do
    tsv="${QC_CHRPLOTS}/${KEY}_per_chrom.tsv"
    [[ -f "$tsv" ]] && tail -n +2 "$tsv" | awk -v k="$KEY" '{print k"\t"$0}'
  done
} > "$CHR_SUMMARY_TSV"
log "Per-chrom table: $CHR_SUMMARY_TSV"

BG_PATHS=(); BG_LABELS=()
for entry in "${KARYOGRAM_BG_LIST[@]}"; do
  BG_LABELS+=("${entry%%:*}")
  BG_PATHS+=("${entry##*:}")
done

if [[ -n "$KARYOGRAM_PY" && ${#BG_PATHS[@]} -gt 0 ]]; then
  FIRST_KEY="${SAMPLE_KEYS[0]}"
  FIRST_GENOME="${SAMPLE_GENOME[$FIRST_KEY]:-hg38}"
  [[ "${FIRST_GENOME,,}" == "hg38" ]] && CS="${CHROM_SIZES_HUMAN:-}" \
                                       || CS="${CHROM_SIZES_MOUSE:-}"
  if [[ -n "${CS:-}" && -f "${CS:-}" ]]; then
    log "  Generating karyogram plots (${#BG_PATHS[@]} samples)..."
    python3 "$KARYOGRAM_PY" \
      --bedgraph "${BG_PATHS[@]}" \
      --labels "${BG_LABELS[@]}" \
      --genome "$FIRST_GENOME" \
      --chrom-sizes "$CS" \
      --outdir "$QC_CHRPLOTS" \
      >> "$MAIN_LOG" 2>&1 \
      && log "  Karyogram plots complete → $QC_CHRPLOTS/" \
      || warn "  Karyogram plots FAILED (check log)"
  else
    warn "CHROM_SIZES not set or missing — skipping karyogram plots"
  fi
else
  warn "Skipping karyogram: no bedGraphs or plot_chrom_coverage.py missing"
fi

log "Phase 2 complete."

# =============================================================================
# Phase 3: deepTools genome-wide QC
# [A] multiBamSummary bins    } parallel background
# [B] multiBamSummary peaks   } parallel background
# [C] plotFingerprint           serial foreground (uses THREADS while A+B run)
# wait A → [D+E] plotCorrelation pearson+spearman parallel
#         → [F] plotPCA bins
# wait B → [G] plotCorrelation peaks
#         → [H] plotPCA peaks
# =============================================================================
log "=== Phase 3: deepTools genome-wide QC ==="

BAM_LIST=("${SAMPLE_BAMS[@]}")
LABELS=("${SAMPLE_KEYS[@]}")

BINS_NPZ="${QC_DT}/multiBamSummary_bins.npz"
PEAKS_NPZ="${QC_DT}/multiBamSummary_peaks.npz"

log "  [A] multiBamSummary bins (background)..."
(
  multiBamSummary bins \
    -b "${BAM_LIST[@]}" \
    --labels "${LABELS[@]}" \
    -p "$THREADS" \
    -o "$BINS_NPZ" \
    --outRawCounts "${QC_MATRICES}/multiBamSummary_bins.tab" \
    >> "$MAIN_LOG" 2>&1 \
    && log "  [A] multiBamSummary (bins) complete" \
    || warn "  [A] multiBamSummary (bins) failed"
) &
BINS_PID=$!

if [[ -f "${CONSENSUS_PEAK:-}" && -s "${CONSENSUS_PEAK:-}" ]]; then
  log "  [B] multiBamSummary peaks (background)..."
  (
    multiBamSummary BED-file \
      --BED "$CONSENSUS_PEAK" \
      -b "${BAM_LIST[@]}" \
      --labels "${LABELS[@]}" \
      -p "$THREADS" \
      -o "$PEAKS_NPZ" \
      --outRawCounts "${QC_MATRICES}/multiBamSummary_peaks.tab" \
      >> "$MAIN_LOG" 2>&1 \
      && log "  [B] multiBamSummary (peaks) complete" \
      || warn "  [B] multiBamSummary (peaks) failed"
  ) &
  PEAKS_PID=$!
else
  PEAKS_PID=""
  warn "  [B] No consensus peaks — skipping multiBamSummary BED-file"
fi

log "  [C] plotFingerprint (serial, foreground)..."
plotFingerprint \
  -b "${BAM_LIST[@]}" \
  --labels "${LABELS[@]}" \
  -p "$THREADS" \
  --minMappingQuality 20 \
  --skipZeros \
  --plotFile "${QC_PLOTS}/fingerprint.png" \
  --outRawCounts "${QC_MATRICES}/fingerprint.tab" \
  --outQualityMetrics "${QC_TABLES}/fingerprint_metrics.tsv" \
  >> "$MAIN_LOG" 2>&1 \
  && log "  [C] plotFingerprint complete" \
  || warn "  [C] plotFingerprint failed"

wait "$BINS_PID" 2>/dev/null || true
log "  [A] bins NPZ wait done."

if [[ -f "$BINS_NPZ" ]]; then
  log "  [D+E] plotCorrelation pearson + spearman (parallel)..."
  for method in pearson spearman; do
    (
      plotCorrelation \
        -in "$BINS_NPZ" \
        --corMethod "$method" --skipZeros \
        --whatToPlot heatmap --colorMap RdYlBu \
        --plotTitle "Sample ${method} correlation (genome-wide bins)" \
        --plotFile "${QC_PLOTS}/correlation_heatmap_${method}.png" \
        --outFileCorMatrix "${QC_MATRICES}/correlation_matrix_${method}.tsv" \
        >> "$MAIN_LOG" 2>&1 \
        && log "  plotCorrelation ($method) complete" \
        || warn "  plotCorrelation ($method) failed"
    ) &
  done
  wait
  log "  [D+E] plotCorrelation done."

  log "  [F] plotPCA (bins)..."
  plotPCA \
    -in "$BINS_NPZ" \
    --plotTitle "PCA — genome-wide bins" \
    --plotFile "${QC_PLOTS}/pca_bins.png" \
    --outFileNameData "${QC_MATRICES}/pca_bins_data.tsv" \
    >> "$MAIN_LOG" 2>&1 \
    && log "  [F] plotPCA (bins) complete" \
    || warn "  [F] plotPCA (bins) failed"
fi

if [[ -n "$PEAKS_PID" ]]; then
  wait "$PEAKS_PID" 2>/dev/null || true
  log "  [B] peaks NPZ wait done."
fi

if [[ -f "$PEAKS_NPZ" ]]; then
  log "  [G] plotCorrelation (peaks)..."
  plotCorrelation \
    -in "$PEAKS_NPZ" \
    --corMethod pearson --skipZeros \
    --whatToPlot heatmap --colorMap RdYlBu \
    --plotTitle "Sample correlation (consensus peaks)" \
    --plotFile "${QC_PLOTS}/correlation_heatmap_pearson_peaks.png" \
    --outFileCorMatrix "${QC_MATRICES}/correlation_matrix_pearson_peaks.tsv" \
    >> "$MAIN_LOG" 2>&1 \
    && log "  [G] plotCorrelation (peaks) complete" \
    || warn "  [G] plotCorrelation (peaks) failed"

  log "  [H] plotPCA (peaks)..."
  plotPCA \
    -in "$PEAKS_NPZ" \
    --plotTitle "PCA — consensus peaks" \
    --plotFile "${QC_PLOTS}/pca_peaks.png" \
    --outFileNameData "${QC_MATRICES}/pca_peaks_data.tsv" \
    >> "$MAIN_LOG" 2>&1 \
    && log "  [H] plotPCA (peaks) complete" \
    || warn "  [H] plotPCA (peaks) failed"
fi

log "Phase 3 complete."

# =============================================================================
# Phase 4: computeMatrix / plotHeatmap / plotProfile + FRiP over consensus
# SKIP_PROFILE=true skips computeMatrix/plotHeatmap/plotProfile
# =============================================================================
log "=== Phase 4: Signal over consensus peaks + FRiP consensus ==="

BW_EXIST=(); BW_LABELS=()
for i in "${!SAMPLE_KEYS[@]}"; do
  bw="${SAMPLE_BIGWIGS[$i]}"
  [[ -f "$bw" ]] && BW_EXIST+=("$bw") && BW_LABELS+=("${SAMPLE_KEYS[$i]}")
done

if [[ "${SKIP_PROFILE:-false}" == "true" ]]; then
  log "  SKIP_PROFILE=true — skipping computeMatrix/plotHeatmap/plotProfile"
else
  if [[ ${#BW_EXIST[@]} -gt 0 && -f "${CONSENSUS_PEAK:-}" && -s "${CONSENSUS_PEAK:-}" ]]; then
    MATRIX="${QC_MATRICES}/signal_over_peaks.gz"
    log "  computeMatrix (${#BW_EXIST[@]} bigWigs × consensus peaks)..."
    computeMatrix reference-point \
      -S "${BW_EXIST[@]}" \
      -R "$CONSENSUS_PEAK" \
      --referencePoint center \
      -b 2000 -a 2000 \
      --skipZeros \
      --samplesLabel "${BW_LABELS[@]}" \
      -p "$THREADS" \
      -o "$MATRIX" \
      >> "$MAIN_LOG" 2>&1 \
      && log "  computeMatrix complete" \
      || warn "  computeMatrix failed"

    if [[ -f "$MATRIX" ]]; then
      plotHeatmap \
        -m "$MATRIX" \
        --colorMap RdYlBu_r \
        --whatToShow "heatmap and colorbar" \
        --plotTitle "Signal over consensus peaks (±2kb)" \
        --heatmapHeight 15 \
        -o "${QC_PLOTS}/heatmap_signal_over_peaks.png" \
        --outFileSortedRegions "${QC_PEAKS}/heatmap_sorted_regions.bed" \
        >> "$MAIN_LOG" 2>&1 || warn "plotHeatmap failed"

      plotProfile \
        -m "$MATRIX" \
        --plotTitle "Average signal over consensus peaks (±2kb)" \
        --perGroup \
        -o "${QC_PLOTS}/profile_signal_over_peaks.png" \
        --outFileNameData "${QC_MATRICES}/profile_signal_over_peaks.tsv" \
        >> "$MAIN_LOG" 2>&1 || warn "plotProfile failed"
    fi
  fi
fi

# FRiP over consensus peaks — parallel background jobs
if [[ -f "${CONSENSUS_PEAK:-}" && -s "${CONSENSUS_PEAK:-}" ]]; then
  FRIP_CONSENSUS_TSV="${QC_TABLES}/frip_consensus.tsv"
  echo -e "sample_id\ttotal_reads\treads_in_consensus_peaks\tfrip_consensus" > "$FRIP_CONSENSUS_TSV"
  _PIDS=()

  log "  FRiP consensus (parallel, max $QC_PARALLEL_JOBS jobs)..."
  for i in "${!SAMPLE_KEYS[@]}"; do
    KEY="${SAMPLE_KEYS[$i]}"
    BAM="${SAMPLE_BAMS[$i]}"
    TMP_FC="${QC_TMP}/frip_consensus_${KEY}.tmp"
    [[ -f "$TMP_FC" ]] && continue

    _wait_slot "$QC_PARALLEL_JOBS"

    (
      if [[ ! -f "$BAM" ]]; then
        printf '%s\tNA\tNA\tNA\n' "$KEY" > "$TMP_FC"
      else
        TOTAL=$(samtools view -c -F 4 "$BAM" 2>/dev/null || echo 0)
        RIP=$(bedtools intersect -u -a "$BAM" -b "$CONSENSUS_PEAK" -ubam 2>/dev/null \
              | samtools view -c 2>/dev/null || echo 0)
        FRIP_C=$(awk "BEGIN{if($TOTAL>0){printf \"%.4f\",$RIP/$TOTAL}else{print \"NA\"}}")
        printf '%s\t%s\t%s\t%s\n' "$KEY" "$TOTAL" "$RIP" "$FRIP_C" > "$TMP_FC"
      fi
    ) &
    _PIDS+=("$!")
  done
  _drain

  for KEY in "${SAMPLE_KEYS[@]}"; do
    TMP_FC="${QC_TMP}/frip_consensus_${KEY}.tmp"
    [[ -f "$TMP_FC" ]] && cat "$TMP_FC" >> "$FRIP_CONSENSUS_TSV" \
      || echo -e "${KEY}\tNA\tNA\tNA" >> "$FRIP_CONSENSUS_TSV"
  done
  rm -f "${QC_TMP}"/frip_consensus_*.tmp
  log "  FRiP consensus: $FRIP_CONSENSUS_TSV"
fi

rmdir "$QC_TMP" 2>/dev/null || true

# =============================================================================
log ""
log "=== Post-alignment QC complete (v3.1.1) ==="
log "Summary    : $SUMMARY_TSV"
log "Warnings   : $WARNINGS_TSV"
log "Karyograms : $QC_CHRPLOTS/"
log "Plots      : $QC_PLOTS/"
log "Matrices   : $QC_MATRICES/"
log ""
log "Output layout:"
log "  ${OUT_DIR}/"
log "    tables/  — qc_summary.tsv  qc_warnings.tsv  fingerprint_metrics.tsv"
log "               frip_consensus.tsv  per_chromosome_reads.tsv"
log "    plots/"
log "      chromosome_coverage/ — *_karyogram.png (per sample)"
log "                             karyogram_all_samples.png (multi-sample)"
log "                             *_100kb.bedGraph  *_per_chrom.tsv"
log "      fingerprint.png  correlation heatmaps  PCA plots"
log "      heatmap_signal_over_peaks.png  profile_signal_over_peaks.png"
log "    matrices/ — raw count matrices  PCA data  profile TSVs"
log "    logs/     — main log"
log "    peak_sets/ — merged_narrow.bed  merged_broad.bed  consensus_peaks.bed"
log "    deeptools/ — multiBamSummary_bins.npz  multiBamSummary_peaks.npz"
