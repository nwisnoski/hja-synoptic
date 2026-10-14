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
