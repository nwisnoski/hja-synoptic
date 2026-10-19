# Archived files

`repository_tidy_manifest.csv` records original paths, archive paths, reasons,
and checksums. `removed_plankton_workflows_20261010.csv` records the earlier
plankton-workflow removal; Git history retains those scripts and figures.
`cleanup_20261014.csv` records the current figure reorganization and removal of
retired figures, unused mothur code, and obsolete housekeeping machinery. It
gives each original path, destination when moved, SHA-256 checksum, and recovery
Git revision. For removed files, recover the original content with
`git show <git_revision>:<source_path>` after recreating its parent directory.
Raw source inputs, completed inference outputs, and checkpoints remain intact.

The historical sequence-download note
`provenance/download_notes/get_seqs_PRJNA534507.sh` is missing from both locations
and local Git history. Its expected checksum remains in the earlier manifest.
The retired housekeeping tool depended on this missing file; its script and
duplicate plan were removed. Do not invent a replacement for the missing note.

## Side-analysis cleanup, October 19, 2026

[cleanup_20261019.csv](cleanup_20261019.csv) inventories 182 retired files and
nine before-cleanup snapshots, with original paths, reasons, and SHA-256 hashes.
[retired_side_analyses_20261019.tar.gz](retired_side_analyses_20261019.tar.gz)
contains every listed file, including uncommitted work and failed-run records.
The compressed members were checked against their hashes before originals were
removed from the active tree.

The archive contains water-optics and integrated-chemistry comparisons,
standalone sediment condition/diversity analyses, transformation enzyme and
network-context tests, detection-conditioned gradients, and stale root-level
figure copies. Their exclusion follows the manuscript's scope; it does not
change the numerical evidence supporting Figures 1–4. Bulk transformations
remain active in supporting 08 and 09. The formula checks formerly embedded in
the retired context analysis now run directly in 09; counts reproduce exactly and mean fractions
agree within floating-point precision.

To inspect or recover work, extract into a separate directory first:

```sh
mkdir -p /tmp/hja-retired-20261019
tar -xzf archive/retired_side_analyses_20261019.tar.gz -C /tmp/hja-retired-20261019
```

Members retain their original repository-relative paths. Review before copying
anything back: snapshot members describe the pre-cleanup versions of files that
remain active. Main figure scripts, selected PDFs, required molecular/network
results, and cohort inputs were checked for unchanged content.

The active supporting scripts were subsequently renumbered consecutively 01–09.
Archive members and historical manifests retain their original filenames.
