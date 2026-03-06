# ============================================================
# m6A Motif and Amino Acid Composition Analysis
# ------------------------------------------------------------
# This script:
# 1. Counts m6A consensus motifs (RRAC) in DNA sequences
# 2. Computes the proportion of P/D/E/R/K/Q residues in peptides
# 3. Generates a scatter plot comparing motif abundance and
#    amino acid composition
#
# Author: Oktay CAKIL
# ============================================================

# ------------------------------------------------------------
# PART 1: Load Required Packages
# ------------------------------------------------------------

library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# PART 2: Define the m6A Consensus Motif
# ------------------------------------------------------------
# m6A consensus sequence: RRAC
# R corresponds to A or G

motif_m6A <- "[AG][AG]AC"

# ------------------------------------------------------------
# PART 3: Count m6A Motifs in DNA Sequences (INPUT: Supplementary table 2)
# ------------------------------------------------------------

protein_peaks$m6A_motif_count <- str_count(
  protein_peaks$dnaSeq,
  motif_m6A
)

# ------------------------------------------------------------
# PART 4: Calculate Proportion of P/D/E/R/K/Q in Peptides
# ------------------------------------------------------------
# These residues are enriched in low-complexity or disordered regions

protein_peaks$PEDRKQ_prop <- 
  str_count(protein_peaks$pepSeq, "[PDEKRQ]") /
  nchar(protein_peaks$pepSeq)

# ------------------------------------------------------------
# PART 5: Generate Unique Gene Names for Duplicated Entries
# ------------------------------------------------------------
# If multiple peaks correspond to the same gene, we append
# an index (e.g., GeneA(1), GeneA(2), etc.)

protein_peaks <- protein_peaks %>%
  group_by(gene) %>%
  mutate(
    
    gene_count = row_number(),
    gene_total = n(),
    
    gene_unique = ifelse(
      gene_total == 1,
      gene,
      paste0(gene, "(", gene_count, ")")
    )
    
  ) %>%
  ungroup() %>%
  select(-gene_count, -gene_total)

# ------------------------------------------------------------
# PART 6: Scatter Plot
# m6A Motif Count vs P/D/E/R/K/Q Proportion
# ------------------------------------------------------------

p <- ggplot(
  protein_peaks,
  aes(x = PEDRKQ_prop, y = m6A_motif_count)
) +
  
  geom_point(
    alpha = 0.6,
    color = "steelblue"
  ) +
  
  # Label genes with high PEDRKQ enrichment
  geom_text_repel(
    data = subset(protein_peaks, PEDRKQ_prop > 0.4),
    aes(label = gene_unique),
    size = 3
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Proportion of P/D/E/R/K/Q residues",
    y = "Number of m6A motifs per peak",
    title = "Relationship between m6A Motif Occurrence and PEDRKQ Enrichment"
  )

print(p)

# ------------------------------------------------------------
# PART 7: Save Output Figures
# ------------------------------------------------------------

# PNG (high resolution for presentations)
ggsave(
  "results/graph_PEDRKQ_vs_m6A.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# SVG (vector format for publications)
ggsave(
  "results/graph_PEDRKQ_vs_m6A.svg",
  plot = p,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# PART 8: Export Annotated Data
# ------------------------------------------------------------

writexl::write_xlsx(
  protein_peaks,
  "results/protein_peaks_m6A.xlsx"
)

