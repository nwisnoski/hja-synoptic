# Archived files

`repository_tidy_manifest.csv` records original paths, archive paths, reasons,
and checksums. `removed_plankton_workflows_20261010.csv` records the earlier
plankton-workflow removal; Git history retains those scripts and figures.
Raw source inputs and completed sequence outputs remain intact.

Preview the housekeeping manifest with `Rscript archive/repository_housekeeping.R`;
`--restore` reverses its recorded moves. The historical sequence-download note
`provenance/download_notes/get_seqs_PRJNA534507.sh` is missing from both locations
and local Git history. Its expected checksum remains recorded; dry runs flag it
and apply/restore stops. Do not invent a replacement.
