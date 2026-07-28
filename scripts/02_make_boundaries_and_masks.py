#!/usr/bin/env python3
"""
Build Zambia province boundaries and the gridded province masks.

This is the only step in the pipeline that needs the geospatial stack (geopandas /
rasterio). Everything downstream consumes its committed outputs - a small GeoJSON and
two small NetCDF masks - so the analysis notebook has no geopandas dependency at all.
That matters for reproducibility: geopandas/GDAL is the most common install failure for
someone cloning this repo, and they should not need it just to re-run the analysis.

Outputs
  boundaries/zambia_provinces_wgs84.geojson   simplified outlines, for plotting only
  boundaries/province_codes.csv               integer code -> province name
  data/processed/zambia_province_mask.nc          province code on the CHIRPS grid
  data/processed/zambia_province_mask_era5grid.nc province code on the ERA5 grid

Note the two masks are NOT interchangeable: CHIRPS cell centres sit on .125 offsets
(21.625, 21.875, ...) while ERA5 sits on .0/.25 (21.5, 21.75, ...). Using one grid's mask
against the other silently shifts every cell by half a cell.
"""
import os
import geopandas as gpd
import numpy as np
import netCDF4 as nc
from rasterio.features import rasterize
from rasterio.transform import from_origin

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

SHP = "boundaries/NSDI_Zambia_Districts_2022.shp"
# Tolerance for the plotting geometry, in degrees. ~0.005 deg is roughly 500 m - far
# below the 0.25 deg (~28 km) analysis grid, so it is invisible on every map here, but
# it takes the province file from ~10 MB to ~100 KB and keeps the repo cloneable.
SIMPLIFY_DEG = 0.005


def main():
    gdf = gpd.read_file(SHP).to_crs(epsg=4326)
    print(f"read {len(gdf)} districts, crs -> EPSG:4326")

    provinces = gdf.dissolve(by="PROVINCE").reset_index()[["PROVINCE", "geometry"]]
    provinces = provinces.sort_values("PROVINCE").reset_index(drop=True)
    codes = list(enumerate(provinces["PROVINCE"], start=1))
    print(f"dissolved to {len(provinces)} provinces")

    with open("boundaries/province_codes.csv", "w") as f:
        f.write("code,province\n")
        for code, name in codes:
            f.write(f"{code},{name}\n")

    # --- masks: rasterise from FULL-precision geometry ------------------------
    # (simplified outlines would move the coastline of each province by up to the
    #  tolerance and could flip which cells count as inside)
    shapes = [(geom, code) for (code, _), geom in zip(codes, provinces.geometry)]

    def build_mask(path, lon, lat_ascending, label):
        res = lon[1] - lon[0]
        transform = from_origin(lon.min() - res / 2, lat_ascending.max() + res / 2,
                                res, res)
        # rasterize fills north-to-south; flip back to match the ascending lat array
        mask = np.flipud(rasterize(shapes,
                                   out_shape=(len(lat_ascending), len(lon)),
                                   transform=transform, fill=0, dtype="int16"))

        ds = nc.Dataset(path, "w", format="NETCDF4_CLASSIC")
        ds.createDimension("lat", len(lat_ascending))
        ds.createDimension("lon", len(lon))
        v = ds.createVariable("lat", "f8", ("lat",)); v[:] = lat_ascending
        v.units = "degrees_north"
        v = ds.createVariable("lon", "f8", ("lon",)); v[:] = lon
        v.units = "degrees_east"
        m = ds.createVariable("province_code", "i2", ("lat", "lon"), zlib=True)
        m[:, :] = mask
        m.long_name = "Zambia province code (0 = outside Zambia)"
        m.province_codes = "; ".join(f"{c}={n}" for c, n in codes)
        ds.close()
        present = sorted(set(np.unique(mask)) - {0})
        print(f"  {label}: {len(lat_ascending)}x{len(lon)}, "
              f"{int((mask >= 1).sum())} cells inside, {len(present)}/10 provinces present")
        assert len(present) == 10, f"{label}: a province vanished from the mask"

    # ERA5 grid: centres at .0/.25
    build_mask("data/processed/zambia_province_mask_era5grid.nc",
               21.5 + 0.25 * np.arange(52),
               np.sort(-7.75 - 0.25 * np.arange(45)), "ERA5 grid")

    # --- simplified outlines for plotting -------------------------------------
    simple = provinces.copy()
    simple["geometry"] = simple.geometry.simplify(SIMPLIFY_DEG, preserve_topology=True)
    out = "boundaries/zambia_provinces_wgs84.geojson"
    if os.path.exists(out):
        os.remove(out)
    simple.to_file(out, driver="GeoJSON")
    print(f"wrote {out} ({os.path.getsize(out)/1e6:.2f} MB, "
          f"simplify tolerance {SIMPLIFY_DEG} deg)")


if __name__ == "__main__":
    main()
