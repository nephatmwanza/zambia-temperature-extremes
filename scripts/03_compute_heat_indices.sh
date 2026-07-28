#!/usr/bin/env bash
#
# Zambia temperature-extreme indices from ERA5 daily maximum temperature.
#
# Companion to the rainfall analysis at
#   https://github.com/nephatmwanza/zambia-rainfall-extremes
# The two deliberately share a period (1981-2025), a baseline (1991-2020), a domain and a
# statistical pipeline, so that "rainfall shows no detectable trend but temperature does"
# is a like-for-like comparison rather than an artefact of different choices.
#
# PERIOD
#   Annual indices: calendar years 1981-2025. 2026 is dropped because the ERA5 record ends
#   in May 2026 and a partial year would bias every annual statistic downward for TXx and
#   upward for nothing in particular - either way it is not comparable.
#   Seasonal indices: ONDJFM seasons 1981/82 - 2025/26, matching the rainfall analysis.
#
# BASELINE AND THE CALENDAR-DAY PERCENTILE
#   TX90p and WSDI are defined against the 90th percentile of the base period computed
#   SEPARATELY FOR EACH CALENDAR DAY, using a 5-day window centred on that day (ETCCDI).
#   This matters enormously here: Zambia's mean daily maximum runs about 25 C in July and
#   31 C in October, so a single fixed threshold would flag most of the hot season and
#   almost none of the cool season. `ydrunpctl,90,5` builds the day-of-year varying
#   threshold; verified values are 27.1 C on 1 July against 33.9 C on 1 October.
#
# WHY WSDI IS ANNUAL ONLY
#   The spell indices in the rainfall repository are computed by slicing a season and
#   re-stamping it onto a synthetic time axis, so that eca_* does not cut the season at
#   31 December. That trick CANNOT be used here: the warm-spell threshold is indexed by
#   day-of-year, so re-stamping would silently compare each day against the wrong day's
#   threshold. WSDI is therefore computed on calendar years only, which is also the ETCCDI
#   definition. The seasonal set is limited to indices that do not depend on runs of
#   consecutive days, where a plain per-season aggregation is unambiguous.
#
# INDICES
#   txm     mean daily maximum temperature                        degC
#   txx     hottest day (maximum daily maximum)                   degC
#   tx90p   share of days above the calendar-day 90th percentile  % of days
#   su35    days with Tmax above 35 C (heat-stress threshold)     days
#   wsdi    warm spell days: >=6 consecutive days above the       days   [annual only]
#           calendar-day 90th percentile
#
#   35 C is used rather than the ETCCDI default "summer day" of 25 C, which would be met on
#   almost every day in Zambia and carry no information. 2.7% of days in the base period
#   exceed 35 C, and it is a threshold with agronomic meaning for maize.
#
set -euo pipefail
cd "$(dirname "$0")/.."

SRC=data/processed/era5_tmax_zambia_1950_present.nc
OUT=data/processed
TMP=data/processed/heat
mkdir -p "$TMP"

FIRST=1981
LAST=2025            # last complete calendar year, and last complete ONDJFM season start
BASE_FIRST=1991
BASE_LAST=2020

run() { cdo -s "$@" 2> >(grep -vE 'HDF5-DIAG|^ +#[0-9]{3}:|major:|minor:|Warning' >&2 || true); }

# CDO refuses to overwrite an existing output file, so clear anything this script owns.
# Without this the script is not re-runnable, which makes it useless for iterating.
rm -f "$OUT"/zambia_{txm,txx,tx90p,su35,wsdi}_annual.nc "$OUT/zambia_heat_annual.nc"
rm -f "$OUT"/zambia_{txm,txx,tx90p,su35}_ondjfm.nc "$OUT/zambia_heat_ondjfm.nc"
rm -f "$TMP"/tmax_c.nc "$TMP/base_c.nc" "$TMP/p90_doy.nc"

echo "=== 1/4  Converting to degC, flipping latitude, trimming the record ==="
# ERA5 stores latitude NORTH-TO-SOUTH (yfirst=-7.75, yinc=-0.25) while CHIRPS and both
# province masks run SOUTH-TO-NORTH. Without invertlat the mask is applied upside down:
# Northern Province statistics would be computed from southern Zambia, silently and with
# entirely plausible-looking numbers. scripts/04 refuses to run on mismatched grids,
# which is how this was caught - but fix it at the source so every product downstream
# shares one convention.
run -f nc4c -z zip -invertlat -subc,273.15 -selyear,${FIRST}/2026 "$SRC" "$TMP/tmax_c.nc"
run -f nc4c -z zip selyear,${BASE_FIRST}/${BASE_LAST} "$TMP/tmax_c.nc" "$TMP/base_c.nc"
run info "$TMP/tmax_c.nc" | sed -n '2p'

echo "=== 2/4  Calendar-day 90th percentile, base ${BASE_FIRST}-${BASE_LAST} (5-day window) ==="
run ydrunmin,5 "$TMP/base_c.nc" "$TMP/_dmin.nc"
run ydrunmax,5 "$TMP/base_c.nc" "$TMP/_dmax.nc"
run ydrunpctl,90,5 "$TMP/base_c.nc" "$TMP/_dmin.nc" "$TMP/_dmax.nc" "$TMP/p90_doy.nc"
rm -f "$TMP/_dmin.nc" "$TMP/_dmax.nc"
echo -n "    threshold on 1 Jul / 1 Oct (Zambia mean, degC): "
echo -n "$(run outputtab,value -fldmean -seltimestep,182 "$TMP/p90_doy.nc" | tail -1) / "
run outputtab,value -fldmean -seltimestep,274 "$TMP/p90_doy.nc" | tail -1

echo "=== 3/4  Annual indices, ${FIRST}-${LAST} ==="
YRS="selyear,${FIRST}/${LAST}"
run -f nc4c -z zip -setname,txm  -yearmean            -$YRS "$TMP/tmax_c.nc" "$OUT/zambia_txm_annual.nc"
run -f nc4c -z zip -setname,txx  -yearmax             -$YRS "$TMP/tmax_c.nc" "$OUT/zambia_txx_annual.nc"
run -f nc4c -z zip -setname,su35 -yearsum -gtc,35     -$YRS "$TMP/tmax_c.nc" "$OUT/zambia_su35_annual.nc"
run -f nc4c -z zip -setname,tx90p -mulc,100 -yearmean -gtc,0 \
    -ydaysub -$YRS "$TMP/tmax_c.nc" "$TMP/p90_doy.nc" "$OUT/zambia_tx90p_annual.nc"
# WSDI must be computed one year at a time. Every eca_* operator reduces its WHOLE input
# to a single timestep - it does not group by calendar year - so handing it the full
# record silently returns one 45-year total (measured: 365.6 days, i.e. the sum, not the
# annual mean). Slicing by real calendar year also keeps day-of-year alignment intact, so
# each day is still compared against its own threshold. eca_hwfi emits the index plus a
# spell-count companion; keep the index.
rm -f "$TMP"/w_*.nc "$TMP/_wsdi.nc"
for y in $(seq $FIRST $LAST); do
  run selyear,$y "$TMP/tmax_c.nc" "$TMP/_y.nc"
  run -setdate,${y}-07-01 -settime,00:00:00 \
      -selname,warm_spell_days_index_wrt_90th_percentile_of_reference_period \
      -eca_hwfi "$TMP/_y.nc" "$TMP/p90_doy.nc" "$TMP/w_${y}.nc"
  rm -f "$TMP/_y.nc"
done
# mergetime takes many inputs, so a prefix operator cannot be chained onto it; rename after
run mergetime "$TMP"/w_*.nc "$TMP/_wsdi.nc"
run -f nc4c -z zip setname,wsdi "$TMP/_wsdi.nc" "$OUT/zambia_wsdi_annual.nc"
rm -f "$TMP"/w_*.nc "$TMP/_wsdi.nc"

run -f nc4c -z zip merge "$OUT"/zambia_{txm,txx,tx90p,su35,wsdi}_annual.nc \
    "$OUT/zambia_heat_annual.nc"

echo "=== 4/4  ONDJFM growing-season indices, ${FIRST}/82-${LAST}/26 ==="
rm -f "$TMP"/s_*.nc
for y in $(seq $FIRST $LAST); do
  y2=$((y+1))
  real="$TMP/real_${y}.nc"
  stamp="setdate,${y}-10-01 -settime,00:00:00"
  run seldate,${y}-10-01,${y2}-03-31 "$TMP/tmax_c.nc" "$real"
  run $stamp -setname,txm   -timmean          "$real" "$TMP/s_txm_${y}.nc"
  run $stamp -setname,txx   -timmax           "$real" "$TMP/s_txx_${y}.nc"
  run $stamp -setname,su35  -timsum -gtc,35   "$real" "$TMP/s_su35_${y}.nc"
  run $stamp -setname,tx90p -mulc,100 -timmean -gtc,0 -ydaysub "$real" "$TMP/p90_doy.nc" \
      "$TMP/s_tx90p_${y}.nc"
  rm -f "$real"
done
for v in txm txx tx90p su35; do
  run -f nc4c -z zip mergetime "$TMP"/s_${v}_*.nc "$OUT/zambia_${v}_ondjfm.nc"
done
rm -f "$TMP"/s_*.nc
run -f nc4c -z zip merge "$OUT"/zambia_{txm,txx,tx90p,su35}_ondjfm.nc \
    "$OUT/zambia_heat_ondjfm.nc"

echo
echo "Done."
ls -lh "$OUT"/zambia_heat_*.nc
echo "annual variables : $(run showname "$OUT/zambia_heat_annual.nc")"
echo "seasonal variables: $(run showname "$OUT/zambia_heat_ondjfm.nc")"
