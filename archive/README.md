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

To preview or reverse the archive operation:

```sh
Rscript analysis/repository_housekeeping.R
Rscript analysis/repository_housekeeping.R --restore
```
