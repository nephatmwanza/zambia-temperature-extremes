#!/usr/bin/env bash
#
# Step 1 - obtain ERA5 daily maximum temperature and clip it to Zambia.
#
# The raw input is a continental daily NetCDF of several GB and is NOT committed. Download
# it into data/raw/, then run this script; everything downstream works on files small
# enough to keep in version control.
#
#   data/raw/era5_tmax_daily_af.nc   ERA5 daily maximum 2 m temperature, 0.25 deg,
#                                    1950-present, African domain, units KELVIN.
#                                    Source: https://cds.climate.copernicus.eu
#                                    (ERA5 daily statistics) or https://climexp.knmi.nl
#
#   Zambia_Administrative_Boundaries_Districts_2020_*.zip
#                                    Zambia district boundaries (NSDI / ZamStats), EPSG:3857.
#
set -euo pipefail
cd "$(dirname "$0")/.."

RAW=data/raw
OUT=data/processed
mkdir -p "$RAW" "$OUT" boundaries

# Zambia spans lon 21.996-33.710, lat -18.077 to -8.272. The box adds a ~0.5 deg buffer so a
# ring of cells survives just outside the border; without it the province mask would clip
# border cells that legitimately overlap Zambia.
BOX="21.5,34.25,-18.75,-7.75"

if [ -f "$RAW/era5_tmax_daily_af.nc" ]; then
  cdo -f nc4c -z zip sellonlatbox,$BOX "$RAW/era5_tmax_daily_af.nc" \
      "$OUT/era5_tmax_zambia_1950_present.nc"
  echo "  ERA5 -> $OUT/era5_tmax_zambia_1950_present.nc"
else
  echo "  SKIP: $RAW/era5_tmax_daily_af.nc not found (see header for source)"
fi

ZIP=$(ls "$RAW"/Zambia_Administrative_Boundaries_Districts_*.zip \
         Zambia_Administrative_Boundaries_Districts_*.zip 2>/dev/null | head -1 || true)
if [ -n "${ZIP:-}" ]; then
  unzip -o "$ZIP" -d boundaries/ >/dev/null
  echo "  unpacked $ZIP -> boundaries/"
else
  echo "  SKIP: no boundary zip found in $RAW/ or repo root"
fi

echo
echo "Done. Next: scripts/02_make_boundaries_and_masks.py"
