# Saved iCAMP runs

| Directory | Status and role | Samples |
| --- | --- | --- |
| `reviewed_pool_castor/` | Completed; primary exploratory view using the reviewed regional pool | 76 soil, porewater, and sediment libraries |
| `sediment_castor/` | Completed; sediment-only regional-pool sensitivity | 34 sediment libraries |
| `catchment_castor/` | Historical; includes excluded libraries and plankton; retained for provenance | Original 100-library four-habitat pool |

Keep canonical `icamp_result.rds` files, manifests/settings, sample identities,
inference diagnostics, and recovery checkpoints. Network-figure tables retain
the site/segment weighting, balancing selections, and descriptive contrasts.
Catchment habitat-affiliation diagnostics apply only to that historical run.

The reviewed-pool result and manifest contain 2,850 assignments and matching
current input hashes. Figure scripts reuse completed results. Distance caches
require matching input hashes, tree-tip order, completion records, and numerical
checks before reuse; checkpoints require the same run identity and versions.

Selected overviews and additional diagnostics are separated in the
[figure guide](../../figures/supplement/icamp/README.md). No inference results
were removed during figure cleanup.
