#!/usr/bin/env bash
# ==========================================================================
# make_trees_PMT.sh   (v2)
# --------------------------------------------------------------------------
# PMT twin of make_trees_XA.sh.
#   * flat trees  -> OUTFILE_DST  (/pnfs .../LED/PMT_outfiles ; analyzer reads there)
#   * QA plots    -> PLOTS_DIR    (local ./Plots under PMT_timing)
# No change to PMTTimingAnalyzer.C is needed; only PMTWaveformTree.h's chain
# pattern must point at OUTFILE_DST (see the one-line edit in the notes).
#
# RUN AFTER the standard SBND setup.
#   ./make_trees_PMT.sh
#   DRYRUN=1 ./make_trees_PMT.sh
# ==========================================================================
set -uo pipefail

# PMT tree-maker: FUNCTION name must equal FILE stem (AnalyzerMakinTree.C -> AnalyzerMakinTree).
MACRO="${MACRO:-/exp/sbnd/data/users/castalyf/PMT_timing/AnalyzerMakinTree.C}"
# shared decoded files (from decode_LED_data.sh; same files for PMT + XA):
DECODED_DST="${DECODED_DST:-/pnfs/sbnd/persistent/users/castalyf/LED/LED_data}"
# where PMT flat trees go -- MUST match the PMTWaveformTree.h chain pattern dir:
OUTFILE_DST="${OUTFILE_DST:-/pnfs/sbnd/persistent/users/castalyf/LED/PMT_outfiles}"
# QA plots (local, fast disk):
PLOTS_DIR="${PLOTS_DIR:-/exp/sbnd/data/users/castalyf/PMT_timing/Plots}"
NPROC="${NPROC:-4}"
DRYRUN="${DRYRUN:-0}"

command -v root >/dev/null || { echo "ERROR: root not found -- source your setup first."; exit 1; }
command -v ifdh >/dev/null || setup ifdhc 2>/dev/null || true
[ -f "$MACRO" ] || { echo "ERROR: PMT tree-maker not found: $MACRO"; exit 1; }
mkdir -p "$PLOTS_DIR"
LOGDIR="$PWD/tree_logs_PMT"; mkdir -p "$LOGDIR"

make_one() {
  local f="$1"
  local base; base="$(basename "$f" .root)"
  local out="$OUTFILE_DST/outfile_ana-${base}.root"
  local log="$LOGDIR/${base}.log"

  if [ -s "$out" ]; then echo "SKIP exists : $base"; return 0; fi
  if [ "$DRYRUN" = 1 ]; then echo "DRYRUN      : $f -> $out"; return 0; fi

  local tmp; tmp="$(mktemp -d)"; cd "$tmp" || return 1
  root -l -b -q "${MACRO}(\"$f\")" > "$log" 2>&1
  local made="outfile_ana-${base}.root"
  if [ -s "$made" ]; then
    if ifdh cp "$made" "$out" 2>>"$log" || cp "$made" "$out" 2>>"$log"; then
      for d in Plots-*; do
        [ -d "$d" ] || continue
        rm -rf "${PLOTS_DIR:?}/$d"; mv "$d" "$PLOTS_DIR/"
      done
      echo "OK          : $base"
    else echo "FAIL copy   : $base  (see $log)"; fi
  else echo "FAIL        : $base  (see $log)"; fi
  cd /; rm -rf "$tmp"
}
export -f make_one
export OUTFILE_DST PLOTS_DIR MACRO LOGDIR DRYRUN

mapfile -t files < <(ls -1 "$DECODED_DST"/*filtered.root 2>/dev/null | grep -v 'outfile_ana' | sort)
echo "Decoded files found: ${#files[@]}   (NPROC=$NPROC)"
[ "${#files[@]}" -eq 0 ] && { echo "None found in $DECODED_DST"; exit 0; }

printf '%s\n' "${files[@]}" | xargs -d '\n' -P "$NPROC" -I{} bash -c 'make_one "$@"' _ {}

echo "----------------------------------------------------------------"
echo "PMT trees -> $OUTFILE_DST"
echo "QA plots  -> $PLOTS_DIR"
echo -n "Total trees: "; ls -1 "$OUTFILE_DST"/outfile_ana-data_EventBuilder*.root 2>/dev/null | wc -l