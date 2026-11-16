#!/usr/bin/env bash
# ==========================================================================
# decode_LED_data.sh
# --------------------------------------------------------------------------
# Decode (+filter) many LED PDS fragment files and stage the outputs on dCache.
# Detector-agnostic: the decode produces the full PDS stream, so the SAME
# decoded files feed BOTH the PMT and the XARAPUCA tree-makers
# (make_trees_PMT.sh / make_trees_XA.sh) -> PMTWaveformAnalyzer.C / XATimingAnalyzer.C.
#
# Each run directory under RAW_SRC holds MANY fragment files, one per
# EventBuilder/art stream:
#     <RAW_SRC>/run_<RUN>/data_EventBuilder*_run<RUN>_*/<RAW_NAME>
# We decode each of those and copy the filtered result to DECODED_DST under a
# unique name that ends in "_filtered.root" (so the tree-maker globs match).
#
# RUN THIS *AFTER* the standard SBND setup, from your larsoft area:
#   setup_container ; sbnd_login ; cd <larsoft> ; source localProducts*/setup ; mrbslp
# then:
#   ./decode_LED_data.sh                 # decode everything
#   MAXFRAG=3 ./decode_LED_data.sh       # quick test: 3 fragments/run
#   DRYRUN=1 ./decode_LED_data.sh        # just print what would run
# Env vars below can all be overridden on the command line.
# ==========================================================================
set -uo pipefail

# ------------------------------- config -----------------------------------
FCL="${FCL:-run_decoders_job.fcl}"            # decode+filter fcl (on FHICL_FILE_PATH)
RAW_SRC="${RAW_SRC:-/pnfs/sbnd/persistent/users/jromeo/jromeo/XA_data/LED}"
RAW_NAME="${RAW_NAME:-xarapucadecoder-art.root}"   # raw fragment filename to decode
# Shared decoded output -- written here, read by BOTH make_trees_PMT.sh and make_trees_XA.sh:
DECODED_DST="${DECODED_DST:-/pnfs/sbnd/persistent/users/castalyf/LED/LED_data}"

RUNS=(${RUNS:-19403 19410 19417 19424 19429 19434 19439 19449 19454 19463})

MAXFRAG="${MAXFRAG:-0}"     # 0 = ALL fragments/run ; N = cap at N fragments/run
NPROC="${NPROC:-4}"         # parallel lar jobs (shared gpvm: keep modest, 4-8)
DRYRUN="${DRYRUN:-0}"       # 1 = list only, do nothing
# --------------------------------------------------------------------------

command -v lar  >/dev/null || { echo "ERROR: 'lar' not found -- source your larsoft setup first."; exit 1; }
command -v ifdh >/dev/null || { echo "NOTE: 'ifdh' not found, trying 'setup ifdhc'..."; setup ifdhc 2>/dev/null || true; }

mkdir -p "$DECODED_DST" 2>/dev/null || true
LOGDIR="$PWD/decode_logs"; mkdir -p "$LOGDIR"

# ---- build work list: one TAB-separated line per fragment ----------------
worklist="$(mktemp)"
for R in "${RUNS[@]}"; do
  rundir="$RAW_SRC/run_${R}"
  if [ ! -d "$rundir" ]; then echo "WARN: missing $rundir" >&2; continue; fi
  mapfile -t frags < <(find "$rundir" -name "$RAW_NAME" 2>/dev/null | sort)
  if [ "$MAXFRAG" -gt 0 ] 2>/dev/null; then frags=("${frags[@]:0:$MAXFRAG}"); fi
  for f in "${frags[@]}"; do
    fragdir="$(basename "$(dirname "$f")")"      # unique: contains run + timestamp
    printf '%s\t%s\t%s\n' "$R" "$f" "$fragdir" >> "$worklist"
  done
done

n="$(wc -l < "$worklist")"
echo "Fragments to decode: $n   (runs: ${RUNS[*]} | MAXFRAG=$MAXFRAG | NPROC=$NPROC)"
if [ "$n" -eq 0 ]; then echo "Nothing to do."; rm -f "$worklist"; exit 0; fi

# ---- per-fragment worker (exported for xargs) ----------------------------
decode_one() {
  IFS=$'\t' read -r R RAW fragdir <<< "$1"
  local out="$DECODED_DST/${fragdir}_decoded_filtered.root"
  local log="$LOGDIR/${fragdir}.log"

  if [ -s "$out" ]; then echo "SKIP exists : $fragdir"; return 0; fi
  if [ "$DRYRUN" = 1 ]; then echo "DRYRUN      : run$R  $RAW -> $(basename "$out")"; return 0; fi

  local tmp; tmp="$(mktemp -d)"; cd "$tmp" || return 1
  if lar -c "$FCL" -s "$RAW" -n -1 > "$log" 2>&1; then
    # Locate the filtered output the fcl produced (name may vary between fcls).
    local prod
    prod="$(ls -S *decoded*filtered*.root 2>/dev/null | head -1)"
    [ -z "$prod" ] && prod="$(ls -S *_filtered.root 2>/dev/null | head -1)"
    [ -z "$prod" ] && prod="$(ls -S *.root 2>/dev/null | grep -v -e decoder_hist -e '_hist' | head -1)"
    if [ -n "$prod" ] && [ -s "$prod" ]; then
      if ifdh cp "$prod" "$out" 2>>"$log"; then echo "OK          : $fragdir"
      elif cp "$prod" "$out" 2>>"$log";       then echo "OK (cp)     : $fragdir"
      else echo "FAIL copy   : $fragdir  (see $log)"; fi
    else
      echo "FAIL output : $fragdir  (no filtered .root -- see $log)"
    fi
  else
    echo "FAIL lar    : $fragdir  (see $log)"
  fi
  cd / ; rm -rf "$tmp"
}
export -f decode_one
export DECODED_DST LOGDIR FCL DRYRUN

# ---- run in parallel ------------------------------------------------------
xargs -a "$worklist" -d '\n' -P "$NPROC" -I{} bash -c 'decode_one "$@"' _ {}
rm -f "$worklist"

echo "----------------------------------------------------------------"
echo "Done. Decoded files staged in: $DECODED_DST"
echo -n "Total decoded (_filtered.root): "
ls -1 "$DECODED_DST"/*_filtered.root 2>/dev/null | wc -l
echo "Per-fragment logs: $LOGDIR/"