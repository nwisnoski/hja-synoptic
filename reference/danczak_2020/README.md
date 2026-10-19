# Danczak et al. transformation reference

Unmodified files from Robert E. Danczak's
[Meta-Metabolome_Ecology repository](https://github.com/danczakre/Meta-Metabolome_Ecology),
associated with [Danczak et al. (2020), *Using metacommunity ecology to understand
environmental metabolomes*](https://doi.org/10.1038/s41467-020-19989-y).
The original script credits Robert Danczak, PNNL.

Retrieved October 18, 2026, at commit
`434e28bd471c83e1eb616f4203d6ada584520d2e` (November 2, 2020):

| File | SHA-256 |
| --- | --- |
| [Transformation_Analysis_by_Sample.R](https://github.com/danczakre/Meta-Metabolome_Ecology/blob/434e28bd471c83e1eb616f4203d6ada584520d2e/Transformation_Analysis_by_Sample.R) | `ba812026852cad57592b8632fa3fe8121213a6317f7c3433e66fbf2f60afb074` |
| [Transformation_Database_07-2020.csv](https://github.com/danczakre/Meta-Metabolome_Ecology/blob/434e28bd471c83e1eb616f4203d6ada584520d2e/Transformation_Database_07-2020.csv) | `340c3a30a593cd336c115e9e7ad5ec78608c94d6c9654dff86296c33adcb68a5` |

Run `analysis/supporting/08_molecular_transformations.R` from the HJA checkout.
It uses the existing sediment presence/absence filter and measured m/z values.
It calculates differences once across retained peaks, then requires both peaks
to be detected in a sample. The inclusive absolute tolerance is 0.000010 mass
units. As in the upstream script, each unordered pair receives the last matching
database label in original file order. Every alternative assignment remains in
`results/fticr_2016/transformations/peak_pair_assignments.csv.gz`; the
`published_assignment` column identifies the label used in the profiles.

The HJA script checks pair identities and complete label-count profiles against
the unmodified upstream script on one complete sediment sample. It also checks
overlapping labels, inclusive boundaries, and a sample without matches. These
checks are recorded in `numerical_validation.csv`.

Results include the full database-by-site count table, peak-pair assignments,
degrees for all detected peaks (including zeros), site summaries, PCA coordinates
and loadings, and input checksums. Figures describe common labels, Hellinger PCA
of pair-count proportions, and the dependence of match counts on possible pairs.
They contain no environmental hypothesis tests. The existing filter and mass
tolerance condition the results; counts are putative connections, not reaction
rates, structural identifications, or evidence of reaction direction. The
last-label rule can favor formula labels over familiar biochemical names.

For a descriptive overview, run `analysis/supporting/09_transformation_summary.R`
using the completed outputs of 12. It calculates the curated motif formula
checks directly from the molecular inventory and pair assignments. Start with
`results/fticr_2016/transformations/transformation_label_catalogue.csv`: it ranks
all reference labels by their mean within-sample fraction of matched pairs and
retains site coverage, distinct peak pairs, alternative-label matches, and mass
errors. `interpretable_motif_summary.csv` gives readable descriptions and the
formula checks for the declared subset of simple motifs. The two new summary
figures show the full label-rank distribution and distinct pairs by assignment
check. These checks are evidence about assignments, not calibrated probabilities
that a biochemical reaction occurred. In the single-candidate category, both
endpoints have one reported formula candidate; the multiple-candidate category
has at least one ambiguous endpoint.
