"""Estimate soil membership in aquatic-site drainage areas from a public DEM.

Run once before 11_soil_stream_localization.R. Requires pysheds, rasterio,
numpy, pandas, and pyproj. No sequencing data are analyzed in this script.
"""
from pathlib import Path
import hashlib
import json
import time

import numpy as np
import pandas as pd
import rasterio
from rasterio.features import shapes
from pyproj import Transformer
from pysheds.grid import Grid
import pysheds

# 1. Settings ------------------------------------------------------------------
root = Path(__file__).resolve().parents[2]
input_dir = root / "data/spatial/soil_routing"
output_dir = root / "data/derived/soil_drainage_2016"
results_dir = root / "results/diversity_2016/soil_stream_localization"
snap_radii_m = [15, 30, 50]
primary_snap_radius_m = 30
soil_offsets_m = [-20, 0, 20]
acceptable_area_ratio = [0.5, 2.0]
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
output_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)
start_time = time.monotonic()

# 2. Condition the DEM with established pysheds routines ------------------------
dem_path = input_dir / "usgs_3dep_10m_dem.tif"
metadata_path = root / "data/derived/microbial_diversity_2016/sample_metadata.csv"
grid = Grid.from_raster(str(dem_path))
dem = grid.read_raster(str(dem_path))
assert np.isfinite(dem).all()
assert np.min(dem) > 0
with rasterio.open(dem_path) as source:
    assert source.crs.to_epsg() == 26910
    assert source.res == (10.0, 10.0)
    transform = source.transform
print("Conditioning 10-m DEM", flush=True)
pit_filled = grid.fill_pits(dem)
depression_filled = grid.fill_depressions(pit_filled)
inflated = grid.resolve_flats(depression_filled)
fdir = grid.flowdir(inflated, dirmap=dirmap, routing="d8")
acc = grid.accumulation(fdir, dirmap=dirmap, routing="d8")
conditioning = pd.DataFrame([{
    "cells": dem.size,
    "depression_modified_cells": int(np.sum(depression_filled > dem)),
    "maximum_fill_m": float(np.max(depression_filled - dem)),
    "mean_fill_m": float(np.mean(depression_filled - dem)),
    "pixel_area_m2": 100,
}])
conditioning.to_csv(results_dir / "dem_conditioning_audit.csv", index=False)

# 3. Preserve all soils, and transform only coordinates with known identities ----
metadata = pd.read_csv(metadata_path, dtype={"site_code": str})
metadata = metadata.loc[metadata.included_10k].copy()
soils = metadata.loc[metadata.habitat == "soil", [
    "sample_id", "site_code", "mapping_status", "soil_latitude", "soil_longitude"
]].copy()
soils["gps_available"] = soils.soil_latitude.notna() & soils.soil_longitude.notna()
soils["spatial_included"] = soils.gps_available & (soils.mapping_status == "MATCHED")
project = Transformer.from_crs(4326, 26910, always_xy=True)
soils["soil_utm_x_m"] = np.nan
soils["soil_utm_y_m"] = np.nan
mapped = soils.spatial_included
xs, ys = project.transform(soils.loc[mapped, "soil_longitude"].to_numpy(),
                           soils.loc[mapped, "soil_latitude"].to_numpy())
soils.loc[mapped, "soil_utm_x_m"] = xs
soils.loc[mapped, "soil_utm_y_m"] = ys
soils["exclusion_reason"] = ""
soils.loc[~soils.gps_available, "exclusion_reason"] = "Missing GPS coordinates"
soils.loc[soils.mapping_status != "MATCHED", "exclusion_reason"] = "Unresolved soil identity"
soils.to_csv(output_dir / "soil_coordinate_audit.csv", index=False)
spatial_soils = soils.loc[mapped].copy()
sites = metadata.loc[metadata.habitat != "soil", [
    "site_code", "site_utm_x_m", "site_utm_y_m", "site_drainage_area_ha",
    "site_segment_number", "site_stream_order"
]].drop_duplicates("site_code").reset_index(drop=True)
assert sites.notna().all().all()

# 4. Delineate each aquatic site's drainage area under three snapping radii ------
# Select the highest-accumulation pixel within each radius, a standard outlet
# snapping rule. Recorded drainage areas provide independent validation.
site_rows = []
membership_rows = []
features = []
height, width = dem.shape
for site_index, site in sites.iterrows():
    print(f"Delineating {site.site_code} ({site_index + 1}/{len(sites)})", flush=True)
    site_x = float(site.site_utm_x_m)
    site_y = float(site.site_utm_y_m)
    center_col, center_row = (~transform) * (site_x, site_y)
    center_col = int(np.floor(center_col))
    center_row = int(np.floor(center_row))
    for radius in snap_radii_m:
        window = int(np.ceil(radius / 10)) + 1
        rows, cols = np.meshgrid(
            np.arange(max(0, center_row - window), min(height, center_row + window + 1)),
            np.arange(max(0, center_col - window), min(width, center_col + window + 1)),
            indexing="ij")
        pixel_x, pixel_y = transform * (cols + 0.5, rows + 0.5)
        distance = np.hypot(pixel_x - site_x, pixel_y - site_y)
        valid = distance <= radius
        candidates = np.where(valid, np.asarray(acc)[rows, cols], -np.inf)
        selected = np.unravel_index(np.argmax(candidates), candidates.shape)
        outlet_row = int(rows[selected])
        outlet_col = int(cols[selected])
        catch = grid.catchment(x=outlet_col, y=outlet_row, fdir=fdir,
                               dirmap=dirmap, xytype="index", routing="d8")
        catch_array = np.asarray(catch, dtype=bool)
        area_ha = float(catch_array.sum() / 100)
        reference_ha = float(site.site_drainage_area_ha)
        ratio = area_ha / reference_ha
        touches_edge = bool(catch_array[0, :].any() or catch_array[-1, :].any()
                            or catch_array[:, 0].any() or catch_array[:, -1].any())
        area_ok = acceptable_area_ratio[0] <= ratio <= acceptable_area_ratio[1]
        row = dict(site)
        row.update(snap_radius_m=radius, snap_distance_m=float(distance[selected]),
                   outlet_row=outlet_row, outlet_col=outlet_col,
                   outlet_x_m=float(pixel_x[selected]), outlet_y_m=float(pixel_y[selected]),
                   dem_drainage_area_ha=area_ha, drainage_area_ratio=ratio,
                   acceptable_area_ratio=area_ok, catchment_touches_dem_edge=touches_edge,
                   routing_qc_pass=area_ok and not touches_edge)
        site_rows.append(row)
        for _, soil in spatial_soils.iterrows():
            for dx in soil_offsets_m:
                for dy in soil_offsets_m:
                    col, rr = (~transform) * (soil.soil_utm_x_m + dx, soil.soil_utm_y_m + dy)
                    col = int(np.floor(col))
                    rr = int(np.floor(rr))
                    assert 0 <= rr < height and 0 <= col < width
                    membership_rows.append({
                        "soil_sample_id": soil.sample_id, "soil_site_code": soil.site_code,
                        "aquatic_site_code": site.site_code, "snap_radius_m": radius,
                        "soil_offset_x_m": dx, "soil_offset_y_m": dy,
                        "contributing": bool(catch_array[rr, col]),
                        "routing_qc_pass": row["routing_qc_pass"],
                    })
        # Keep drainage polygons for review; R simplifies them only for plotting.
        if radius == primary_snap_radius_m:
            for geometry, value in shapes(catch_array.astype("uint8"),
                                           mask=catch_array, transform=transform):
                features.append({"type": "Feature", "geometry": geometry,
                                 "properties": {"site_code": site.site_code,
                                                "dem_area_ha": area_ha,
                                                "routing_qc_pass": row["routing_qc_pass"]}})

pd.DataFrame(site_rows).to_csv(output_dir / "aquatic_catchment_audit.csv", index=False)
membership = pd.DataFrame(membership_rows)
membership.to_csv(output_dir / "soil_site_membership_sensitivity.csv", index=False)
primary = membership.loc[(membership.snap_radius_m == primary_snap_radius_m)
                        & (membership.soil_offset_x_m == 0)
                        & (membership.soil_offset_y_m == 0)].copy()
robust = membership.groupby(["soil_sample_id", "aquatic_site_code"]).agg(
    contributing_fraction=("contributing", "mean"),
    all_routing_qc_pass=("routing_qc_pass", "all"))
primary = primary.merge(robust, on=["soil_sample_id", "aquatic_site_code"], validate="one_to_one")
primary["membership_stable"] = primary.contributing_fraction.isin([0, 1])
primary["robust_spatial_included"] = primary.membership_stable & primary.all_routing_qc_pass
primary.to_csv(output_dir / "soil_site_membership.csv", index=False)
with open(output_dir / "aquatic_catchments_30m.geojson", "w") as handle:
    json.dump({"type": "FeatureCollection", "crs": {
        "type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::26910"}},
        "features": features}, handle)
checksums = []
for path in [dem_path, metadata_path]:
    checksums.append({"path": str(path.relative_to(root)),
                      "md5": hashlib.md5(path.read_bytes()).hexdigest()})
pd.DataFrame(checksums).to_csv(results_dir / "routing_input_checksums.csv", index=False)
with open(results_dir / "routing_session.json", "w") as handle:
    json.dump({"pysheds_version": pysheds.__version__, "elapsed_seconds": time.monotonic() - start_time,
               "routing": "D8", "snap_radii_m": snap_radii_m,
               "soil_offsets_m": soil_offsets_m, "acceptable_area_ratio": acceptable_area_ratio},
              handle, indent=2)
print(primary.groupby("soil_sample_id").contributing.sum().to_string(), flush=True)
print("Routing complete", flush=True)
