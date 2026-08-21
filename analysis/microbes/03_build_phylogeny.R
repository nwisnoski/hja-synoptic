# Build a master phylogeny for every screened bacterial and archaeal ASV.
# Run 01_prepare_diversity.R first, then run this file from the repository root.

library(Biostrings)
library(DECIPHER)
library(ape)

fasttree <- Sys.getenv(
  "FASTTREE_BIN",
  unset = "/Users/nawis/opt/FastTree/FastTree"
)
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
asv_file <- "data/derived/microbial_diversity_2016/asv_metadata_screened.csv"
output_dir <- "results/phylogeny_2016"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use the full taxonomically screened ASV catalog, independent of rarefaction.
sequences <- read.csv(asv_file, stringsAsFactors = FALSE)

if (anyNA(sequences$asv_id) || anyDuplicated(sequences$asv_id)) {
  stop("Some ASVs are missing sequences or have duplicate sequence records.")
}
if (any(grepl("[^ACGT]", sequences$sequence))) {
  stop("An ASV sequence contains an unexpected nucleotide character.")
}

dna <- DNAStringSet(sequences$sequence)
names(dna) <- sequences$asv_id

unaligned_file <- file.path(output_dir, "asv_sequences_screened.fasta")
alignment_file <- file.path(output_dir, "asv_alignment_screened.fasta")
tree_file <- file.path(output_dir, "asv_tree_screened_fasttree.nwk")
log_file <- file.path(output_dir, "fasttree_screened.log")

writeXStringSet(dna, unaligned_file, format = "fasta")

# Align the short 16S V4 sequences before estimating the tree. For this
# 74,782-sequence catalog, skip DECIPHER's optional global refinement passes:
# they require a new all-by-all distance matrix after the progressive alignment.
alignment <- AlignSeqs(
  dna,
  iterations = 0,
  refinements = 0,
  processors = threads,
  verbose = TRUE
)
if (length(unique(width(alignment))) != 1L) {
  stop("The multiple sequence alignment does not have a common width.")
}
writeXStringSet(alignment, alignment_file, format = "fasta")

# Use a nucleotide GTR model, gamma-rescaled branch lengths, and the faster
# neighbor-joining search appropriate for a data set of more than 50,000 tips.
if (!file.exists(fasttree)) stop("FastTree was not found at: ", fasttree)
Sys.setenv(OMP_NUM_THREADS = threads)
status <- system2(
  fasttree,
  args = c(
    "-nt", "-gtr", "-gamma", "-fastest", "-pseudo",
    "-seed", "2016",
    "-log", log_file,
    "-out", tree_file,
    alignment_file
  )
)
if (status != 0L || !file.exists(tree_file)) {
  stop("FastTree failed; inspect: ", log_file)
}

# Confirm that the Newick tree contains every screened ASV exactly once.
tree <- read.tree(tree_file)
if (anyDuplicated(tree$tip.label) ||
    !setequal(tree$tip.label, sequences$asv_id)) {
  stop("The FastTree tip labels do not match the screened ASV catalog.")
}

tree_summary <- data.frame(
  asvs = length(dna),
  minimum_sequence_length = min(width(dna)),
  median_sequence_length = median(width(dna)),
  maximum_sequence_length = max(width(dna)),
  alignment_columns = unique(width(alignment)),
  tree_tips = length(tree$tip.label),
  fasttree_model = "GTR+CAT followed by Gamma20 branch-length rescaling",
  rooted = FALSE
)
write.csv(
  tree_summary,
  file.path(output_dir, "phylogeny_summary.csv"),
  row.names = FALSE
)

message("Phylogeny complete: ", tree_file)
