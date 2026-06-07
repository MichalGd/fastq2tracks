#!/bin/bash
# fastq2tracks v3.0.2 — MACS2 peak calling — always runs BOTH narrow AND broad
# Usage: bash scripts/macs2_peaks.sh <ip.bam> <ctrl.bam|none> <outDir> <mode> <genome_key> [sample_name]
set -euo pipefail
_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.sh not found. Export F2T_CONFIG or pass --config to fastq2tracks.sh." >&2
            exit 1
        }
    fi
}
_load_config

IP_BAM="$1"; CTRL_BAM="$2"; OUT_DIR="$3"
MODE="${4:-both}"; GENOME_KEY="${5:-hg38}"; SAMPLE="${6:-$(basename "$IP_BAM" .bam)}"
mkdir -p "${OUT_DIR}/narrow" "${OUT_DIR}/broad"
[[ "$GENOME_KEY" == "hg38" ]] && GSIZE="$MACS2_GENOME_HG38" || GSIZE="$MACS2_GENOME_MM39"
CTRL_FLAG=""; [[ "$CTRL_BAM" != "none" && -f "$CTRL_BAM" ]] && CTRL_FLAG="-c $CTRL_BAM"

run_macs2() {
    local peak_type="$1" sub="$2" extra_flags="$3"
    echo "[MACS2] $SAMPLE mode=$peak_type"
    macs3 callpeak \
        -t "$IP_BAM" $CTRL_FLAG \
        -f BAM -g "$GSIZE" -n "$SAMPLE" \
        --outdir "${OUT_DIR}/${sub}" \
        -q "${MACS2_QVALUE}" $extra_flags \
        2>"${OUT_DIR}/${sub}/${SAMPLE}_macs2_${peak_type}.log"
}
if [[ "${MODE,,}" == "none" ]]; then
    echo "MACS2 skipped for $SAMPLE (mode=none)"; exit 0
fi
run_macs2 "narrowPeak" "narrow" ""
run_macs2 "broadPeak"  "broad"  "--broad --broad-cutoff ${MACS2_BROAD_CUTOFF}"
