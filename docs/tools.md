# Tools in qpgentools

`qpgentools` is a collection of collection of C++ tools that are designed to help process modified version of PLINK2 files to rapidly perform quality control, normalization, and evaluation of analysis pipelines.
Currently, the following tools are available

* xQTL Analysis Tools
  * [`qpgentools pair-assoc`](tools/pair_assoc.md): Perform rapid association tests between selected variant-trait pairs.
  * [`qpgentools pair-prs`](tools/pair_prs.md): Perform rapid calculation of polygenic risk scores (PRS) based on selected variant-trait pairs
  * [`qpgentools rect-assoc`](tools/rect_assoc.md): Perform rapid association tests between all pairs between selected subset of variants and selected subset of traits 
  * [`qpgentools match-prs-pheno](tools/match_prs_pheno.md): Perform sample identity matching based on PRS and phenotype data.
* Sequence Data Tools
  * [`qpgentools match-plp-geno`](tools/match_plp_geno.md): Identify the best matching samples based on sequence pileups
  * [`qpgentools join-plp-geno`](tools/join_plp_geno.md): Join the sequence pileups and genotypes for visual inspection
* Genotype Extraction Tools
  * [`qpgentools geno2tsv`](tools/geno2tsv.md): Extract genotypes from PLINK2 files to TSV format.

!!! note
    This documentation contains only a subset of primary tools implemented `qpgentools`.
    Other tools not documented here are considered non-primary tools that have limited support, and they are NOT supported officially.