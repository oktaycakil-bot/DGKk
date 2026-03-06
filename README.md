# DGKk
Cakil et al.

# Sequence Motif Analysis

R scripts for exploratory analysis of sequence features in peptide and DNA datasets.

The repository includes:
- **tripeptide_distribution.R** – detects enriched unordered tripeptides in protein sequences using sliding windows and visualizes their hydrophobicity, frequency, and charge.
- **m6A_count_protein_peaks.R** – quantifies the m6A consensus motif (RRAC) in DNA sequences and compares it to the proportion of P/D/E/R/K/Q residues in corresponding peptides.

Outputs include publication-ready plots (PNG/SVG) and processed tables.

## Requirements
R packages: `ggplot2`, `dplyr`, `stringr`, `ggrepel`, `writexl`, `scales`.

## Author
Oktay
