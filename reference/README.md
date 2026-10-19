# Reference files

## Molecular transformations

[danczak_2020](danczak_2020/README.md) holds the pinned, unmodified transformation
database and per-sample R script used by `analysis/supporting/08_molecular_transformations.R`.

## Taxonomy

Place DADA2-formatted taxonomy training files here. 

- `silva_nr99_v138.2_toGenus_trainset.fa.gz`
- `silva_v138.2_assignSpecies.fa.gz`

The large FASTA files are intentionally ignored by Git. Taxonomy can also be added
after ASV inference because the exact sequences are preserved in
`asv_sequences.fasta`.
