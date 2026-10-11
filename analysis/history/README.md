# Historical project records

[PROJECT_HISTORY.md](PROJECT_HISTORY.md) preserves the dated October 2-10
continuity records extracted intact from the root AGENTS.md during the
October 10 repository tidy-up. Some records describe missing local files,
earlier sample counts, failed jobs, or plans that were later completed.
Use dates to interpret them; they do not override current instructions.

Current entry points are the [analysis roadmap](../README_ANALYSIS_PLAN.md),
[manuscript outline](../SOIL_POREWATER_SEDIMENT_MANUSCRIPT.md), and
[root instructions](../../AGENTS.md). The
[sediment deployment record](../cluster/README_SEDIMENT_DEPLOYMENT.md) documents
job 432228's approved transfer/submission, not its current queue state.

The older [figure/commit review proposal](../REVIEW_2016_FIGURES_AND_COMMITS.md)
is retained at its original path to preserve its relative figure links. Its
October 4 proposal is historical, not a current staging list or authorization
to commit. The [sediment-only manuscript alternative](../SEDIMENT_ONLY_MANUSCRIPT.md)
likewise remains available, with a pointer to the current comparative outline.

Scientific scripts, data, result files, and figures were left in their existing
locations. The tidy-up corrected current navigation and stale status statements
without changing analyses or broadening ignore rules.

Verification compared SHA-256 hashes for all 724 inventoried scientific/code
files. All 723 analysis/data/output files are unchanged. The sole code edit is
`analysis/repository_housekeeping.R`: dry-run previews now retain missing
provenance as a flagged row and warning, while apply/restore still fails.
The archive's download note is missing and unavailable from local Git history;
the original manifest/checksum are preserved. The other four archive entries
pass their recorded checksums. Local Markdown links resolved, and
`git diff --check` passed. Existing uncommitted analysis and deployment
changes were preserved. No commit, push, or remote command was performed.

## October 10 plankton-workflow removal

`removed_plankton_workflows_20261010.csv` records the 13 obsolete scripts and
47 figures removed at Nathan's request. Each row includes the pre-removal
SHA-256 and recovery commit `908d201`. Use `git show 908d201:path/to/file` to
inspect an original; do not restore obsolete workflows as current manuscript
analyses. Current first-detection outputs have been regenerated on the reviewed
three-habitat cohort. Raw inputs and completed inference were preserved.
