#!/bin/bash
# Runs the validation study: 20 simulated landscapes for each of three
# scenarios, then summarises them.
#   A. valleys, analysed with the true resistance surface and smooth moisture
#   B. the same landscapes, analysed with the slope-only surface and with
#      moisture patches on the land between meadows
#   C. valleys with a stronger warming, analysed as in B
# Usage (from the package root):
#   data-raw/validation/run_validation.sh <output folder> [Rscript] [parallel jobs]
set -e
out=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); R=${2:-Rscript}; jobs=${3:-3}
mkdir -p "$out"
job() {
  scen=$1; seed=$2; d="$out/$scen/rep$(printf %03d "$seed")"; mkdir -p "$out/$scen"
  if [ "$scen" = valleys ]; then
    RESIST=true "$R" data-raw/validation/run_replicate.R "$seed" "$d" "$scen" > "$d.A.log" 2>&1
    REUSE=1 RESIST=slope PATCHES=1 PATCH_SIZE=1500 "$R" data-raw/validation/run_replicate.R "$seed" "$d" "$scen" > "$d.B.log" 2>&1
  else
    RESIST=slope PATCHES=1 PATCH_SIZE=1500 "$R" data-raw/validation/run_replicate.R "$seed" "$d" "$scen" > "$d.C.log" 2>&1
  fi
}
export -f job; export out R
for seed in $(seq 1 20); do echo "valleys $seed"; echo "valleys_warmer $seed"; done |
  xargs -P "$jobs" -L 1 bash -c 'job "$0" "$1"'
"$R" data-raw/validation/summarise.R "$out"
