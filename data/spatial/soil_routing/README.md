# Reference terrain for the 2016 soil drainage analysis

Downloaded October 5, 2026. These reference layers contain no later microbial
survey observations. Keep the frozen terrain extract to reproduce the routing.

`usgs_3dep_10m_dem.tif` is a 2,100 x 1,700 float32 elevation extract from the
[USGS 3DEP bare-earth elevation service](https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer).
Its bounding box is 553000,4890000,574000,4907000 in NAD83 / UTM Zone 10N
(EPSG:26910), with 10-m output pixels. The multi-resolution service's native
source resolution can vary; 10 m describes the exported grid. This is a
current terrain reference, not the original 2008 LiDAR model used by Ward et al.
It is used to estimate drainage membership, not to reconstruct 2016 discharge.
The elevation service reports elevations in meters. The accompanying
`usgs_dem_download_metadata.json` records the returned extent and temporary
download URL; that URL may expire.

To request another extract, use the service's `exportImage` endpoint with:

```text
bbox=553000,4890000,574000,4907000
bboxSR=26910
imageSR=26910
size=2100,1700
format=tiff
pixelType=F32
renderingRule={"rasterFunction":"None"}
f=json
```

A new service extract may differ from the frozen input. Input MD5 hashes are
recorded in `results/diversity_2016/soil_stream_localization/`.

`hja_gaged_watersheds.geojson` contains 12 publicly served gaged watershed
polygons from the Andrews Forest research program's
[studies15_1 layer 76](https://cypress.forestry.oregonstate.edu/arcgis/rest/services/hja/studies15_1/MapServer/76).
The query was `where=1=1`, `outFields=*`, `returnGeometry=true`, `outSR=4326`,
and `f=geojson`. This is a partial watershed inventory, not a full catchment
map for all sampling sites. Three mapped soils fall within Watershed 1.

The local `data/_ReadMe_Geometry.docx` establishes UTM Zone 10N for the
2016 stream centerlines but does not name the horizontal datum. The current
analysis assigns NAD83 / UTM Zone 10N to those coordinates, consistent with
the public Andrews GIS. This datum assumption is documented rather than
treated as verified original metadata. GPS coordinates are interpreted as
WGS84 longitude/latitude; their original datum is also unspecified. The
20-m coordinate-offset sensitivity assesses local assignment instability,
but does not cover every possible mapping error.

Reference: Ward et al. (2019), [Co-located contemporaneous mapping of
morphological, hydrological, chemical, and biological conditions in a 5th-order
mountain stream network](https://doi.org/10.5194/essd-11-1567-2019).
