#!/usr/bin/env bash
# Cache ENA run-level metadata for every candidate CRC cohort in config/cohorts.tsv.
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CFG="$ROOT/config/cohorts.tsv"
OUT="$ROOT/data/metadata"
LOG="$ROOT/logs/01_fetch_metadata.log"
mkdir -p "$OUT" "$ROOT/logs"

FIELDS="run_accession,sample_accession,study_accession,library_strategy,instrument_platform,read_count,base_count,fastq_bytes,sample_title,experiment_title"
: > "$LOG"

tail -n +2 "$CFG" | while IFS=$'\t' read -r study accession repo region design notes; do
  [ -z "${accession:-}" ] && continue
  dest="$OUT/${study}__${accession}.tsv"
  url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=${FIELDS}&format=tsv&limit=0"
  printf '%-22s %-14s ' "$study" "$accession"
  if curl -fsS --max-time 120 "$url" -o "$dest.tmp" && [ -s "$dest.tmp" ]; then
    mv "$dest.tmp" "$dest"
    rows=$(( $(wc -l < "$dest") - 1 ))
    wgs=$(awk -F'\t' 'NR>1 && $4=="WGS"{c++} END{print c+0}' "$dest")
    amp=$(awk -F'\t' 'NR>1 && $4=="AMPLICON"{c++} END{print c+0}' "$dest")
    echo "OK   runs=$rows  WGS=$wgs  AMPLICON=$amp" | tee -a "$LOG"
  else
    rm -f "$dest.tmp"
    echo "FAIL (no data / unresolved accession)" | tee -a "$LOG"
  fi
  sleep 2
done

echo "Wrote metadata to $OUT"
