# Repository archive

Files here are preserved, not deleted. Their original paths, current paths,
reasons, and MD5 checksums are recorded in `repository_tidy_manifest.csv`.

The applied set is deliberately small:

- fuller raw environmental exports formerly named `copy.csv`;
- a superseded one-off sequence-download note;
- historical mothur cluster-output logs.

The legacy mothur R Markdown workflow, its OTU tables, scripts, and figures
remain in their original locations because they still form a connected
historical analysis. Active DADA2 files and intermediates are never housekeeping
targets.

October 10, 2026 checkout audit: the manifest's historical download note
`archive/provenance/download_notes/get_seqs_PRJNA534507.sh` is missing from
both its source and archive paths and is not recoverable from available local
Git history. The other four entries pass their recorded checksum checks.
The manifest and expected checksum are retained; no replacement was invented.
The dry run reports this as a `missing` row and a warning. Apply/restore still
fails if a required file is missing or does not match its checksum.

To preview or reverse the archive operation:

```sh
Rscript analysis/repository_housekeeping.R
Rscript analysis/repository_housekeeping.R --restore
```
