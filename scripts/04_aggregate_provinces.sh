#!/usr/bin/env bash
#
# Aggregate a gridded index file to province and national means.
#
# Uses cdo fldmean, which applies proper cos(latitude) area weighting. Zambia spans
# ~8S-18S, so cell area varies ~5% north to south; a plain unweighted mean would bias
# province and national figures. Do not replace this with a numpy .mean().
#
# The mask MUST be on the same grid as the index file. CHIRPS cell centres sit on .125
# offsets and ERA5 on .0/.25, so the two masks are not interchangeable - using the wrong
# one shifts every cell by half a cell without raising an error. The check below refuses
# to run if the grids differ.
#
# Usage (all arguments optional; defaults do the annual heat aggregation):
#   scripts/04_aggregate_provinces.sh [INDEX_NC] [MASK_NC] [OUT_CSV] [YEAR_COLUMN] [VARS...]
#
# Examples:
#   scripts/04_aggregate_provinces.sh                       # annual heat indices
#   scripts/04_aggregate_provinces.sh \
#       data/processed/zambia_heat_ondjfm.nc \
#       data/processed/zambia_province_mask_era5grid.nc \
#       data/processed/zambia_heat_ondjfm_by_province.csv \
#       season_start_year txm txx tx90p su35
#
set -euo pipefail
cd "$(dirname "$0")/.."

IDX=${1:-data/processed/zambia_heat_annual.nc}
MASK=${2:-data/processed/zambia_province_mask_era5grid.nc}
OUT_CSV=${3:-data/processed/zambia_heat_annual_by_province.csv}
YEAR_COL=${4:-year}
shift 4 2>/dev/null || true
VARS=${*:-"txm txx tx90p su35 wsdi"}

TMP=data/processed/agg_tmp
mkdir -p "$TMP"

run() { cdo -s "$@" 2> >(grep -vE 'HDF5-DIAG|^ +#[0-9]{3}:|major:|minor:|Warning' >&2 || true); }

# --- refuse to aggregate with a mismatched mask -----------------------------------
gi=$(run griddes "$IDX" | grep -E '^(xsize|ysize|xfirst|yfirst|xinc|yinc)' | tr -d ' ')
gm=$(run griddes "$MASK" | grep -E '^(xsize|ysize|xfirst|yfirst|xinc|yinc)' | tr -d ' ')
if [ "$gi" != "$gm" ]; then
  echo "ERROR: index grid does not match mask grid." >&2
  echo "--- index ($IDX) ---" >&2; echo "$gi" >&2
  echo "--- mask  ($MASK) ---" >&2; echo "$gm" >&2
  exit 1
fi

# mask variable is province_code; rename so cdo arithmetic lines up with the index vars
run chname,province_code,m "$MASK" "$TMP/m.nc"

emit_region() {   # $1 = region label, $2 = file holding a 1/0 mask named m
  local label="$1" maskfile="$2"
  for v in $VARS; do
    run -selname,"$v" "$IDX" "$TMP/v.nc"
    # ifthen keeps cells where mask is nonzero/non-missing, then fldmean area-weights them
    run fldmean -ifthen "$maskfile" "$TMP/v.nc" "$TMP/agg.nc"
    run outputtab,year,value "$TMP/agg.nc" \
      | awk -v L="$label" -v V="$v" 'NR>1 && NF>=2 {print L","V","$1","$2}'
    rm -f "$TMP/v.nc" "$TMP/agg.nc"
  done
}

{
  echo "region,index,${YEAR_COL},value"

  # national: every cell inside Zambia (code >= 1)
  run -setname,m -gec,1 "$TMP/m.nc" "$TMP/mask_nat.nc"
  emit_region "Zambia" "$TMP/mask_nat.nc"

  # per province
  while IFS=, read -r code name; do
    [ "$code" = "code" ] && continue
    [ -z "$code" ] && continue
    rm -f "$TMP/mask_p.nc"
    run -setname,m -eqc,"$code" "$TMP/m.nc" "$TMP/mask_p.nc"
    emit_region "$name" "$TMP/mask_p.nc"
  done < boundaries/province_codes.csv
} > "$OUT_CSV"

rm -rf "$TMP"

echo "Wrote $OUT_CSV"
wc -l "$OUT_CSV"
