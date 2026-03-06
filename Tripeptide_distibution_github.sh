# ============================================================
# Tripeptide Enrichment Analysis in Protein Sequences
# ------------------------------------------------------------
# This script identifies enriched unordered tripeptides
# within sliding windows across peptide sequences and
# visualizes their enrichment relative to hydrophobicity.
#
# Author: Oktay CAKIL
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# PART 0: Load Data
# ------------------------------------------------------------

# Read peptide sequence file
# Make sure the separator is ";" and the file contains a header
seq_peptide <- read.csv(
  "data/peptide_sequence_list.csv",
  sep = ";",
  header = TRUE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Utility function: sort amino acids within a motif
# Example: 'PEE' -> 'EEP'
# This allows unordered motif comparison
# ------------------------------------------------------------
sort_motif_aa <- function(motif) {
  aas <- unlist(strsplit(motif, ""))
  paste(sort(aas), collapse = "")
}

# ------------------------------------------------------------
# PART 1: Identify enriched tripeptides in sliding windows
# ------------------------------------------------------------
find_enriched_tripeptides_window <- function(sequences, window_size = 30) {
  
  k <- 3
  all_tripeptides <- list()
  
  # Generate ALL unique unordered tripeptide motifs
  all_unordered_motifs_in_seq <- unique(
    do.call(c, lapply(sequences, function(s) {
      
      n_s <- nchar(s)
      if (n_s < k || is.na(s)) return(NULL)
      
      raw_tripeptides <- substring(s, 1:(n_s - k + 1), k:n_s)
      sapply(raw_tripeptides, sort_motif_aa)
      
    }))
  )
  
  for (seq in sequences) {
    
    if (is.na(seq) || is.null(seq)) next
    
    seq <- as.character(seq)
    n <- nchar(seq)
    if (n < k) next
    
    # Slide a window across the sequence
    for (start_pos in 1:(n - window_size + 1)) {
      
      window_seq <- substring(seq, start_pos, start_pos + window_size - 1)
      
      # Generate all overlapping tripeptides in the window
      n_win <- nchar(window_seq)
      window_raw_tripeptides <- substring(window_seq, 1:(n_win - k + 1), k:n_win)
      
      # Convert to unordered motifs
      window_unordered_tripeptides <- sapply(window_raw_tripeptides, sort_motif_aa)
      
      # Count occurrences of each motif in the window
      for (motif_unordered in all_unordered_motifs_in_seq) {
        
        count_in_window <- sum(window_unordered_tripeptides == motif_unordered)
        
        # Enrichment threshold (>=2 occurrences in window)
        if (count_in_window >= 2) {
          
          all_tripeptides <- append(all_tripeptides, list(list(
            Motif = motif_unordered,
            Count = count_in_window,
            Sequence = seq,
            Start_Window = start_pos
          )))
          
        }
      }
    }
  }
  
  return(all_tripeptides)
}

# ------------------------------------------------------------
# PART 2: Data Preparation
# ------------------------------------------------------------

sequences <- seq_peptide$seq

# Identify enriched tripeptides
enriched_repeats_list <- find_enriched_tripeptides_window(
  sequences,
  window_size = 30
)

if (length(enriched_repeats_list) == 0) {
  stop("No unordered tripeptides detected in 30 AA windows. Consider relaxing enrichment criteria.")
}

# Convert list to dataframe
enriched_repeats_df <- do.call(rbind.data.frame, enriched_repeats_list)
enriched_repeats_df$Count <- as.numeric(enriched_repeats_df$Count)

# ------------------------------------------------------------
# PART 3: Motif Counting and Aggregation
# ------------------------------------------------------------

# Count the number of windows containing each motif
motif_counts <- table(enriched_repeats_df$Motif)

total_windows_with_enrichment <- sum(motif_counts)

motif_stats <- data.frame(
  Motif = names(motif_counts),
  Count_Windows = as.numeric(motif_counts),
  Frequency = as.numeric(motif_counts) / total_windows_with_enrichment
)

# Motif length
motif_stats$Length <- 3

# ------------------------------------------------------------
# PART 3.5: Hydrophobicity Calculation
# ------------------------------------------------------------

# Kyte-Doolittle hydrophobicity scale
aa_hydro <- c(
  "I"=4.5, "V"=4.2, "L"=3.8, "F"=2.8, "C"=2.5,
  "M"=1.9, "A"=1.8, "G"=-0.4, "T"=-0.7, "S"=-0.8,
  "W"=-0.9, "Y"=-1.3, "P"=-1.6, "H"=-3.2,
  "Q"=-3.5, "N"=-3.5, "E"=-3.5, "D"=-3.5,
  "K"=-3.9, "R"=-4.5
)

# Compute average hydrophobicity of a tripeptide
calculate_hydrophobicity <- function(motif) {
  
  if (nchar(motif) != 3) return(NA)
  
  aas <- unlist(strsplit(motif, ""))
  score <- sum(aa_hydro[aas])
  
  return(score / 3)
}

motif_stats$Hydrophobicity <- sapply(
  as.character(motif_stats$Motif),
  calculate_hydrophobicity
)

# ------------------------------------------------------------
# PART 3.6: Net Charge Calculation
# ------------------------------------------------------------

aa_charges <- c(
  "K" = 1,
  "R" = 1,
  "D" = -1,
  "E" = -1
)

calculate_net_charge <- function(motif) {
  
  aas <- unlist(strsplit(motif, ""))
  
  net_charge <- sum(
    sapply(aas, function(aa)
      ifelse(aa %in% names(aa_charges), aa_charges[aa], 0))
  )
  
  return(net_charge)
}

motif_stats$NetCharge <- sapply(
  as.character(motif_stats$Motif),
  calculate_net_charge
)

# ------------------------------------------------------------
# PART 4: Scatter Plot
# Hydrophobicity vs Enrichment Frequency
# ------------------------------------------------------------

p_scatter_hydro <- ggplot(
  motif_stats,
  aes(x = Hydrophobicity, y = Frequency, label = Motif)
) +
  
  geom_point(
    aes(fill = factor(NetCharge)),
    shape = 21,
    colour = "black",
    size = 4,
    alpha = 0.8
  ) +
  
  # Label highly enriched or strongly hydrophobic motifs
  geom_text(
    data = subset(
      motif_stats,
      abs(Hydrophobicity) > 2.5 |
        Frequency > quantile(Frequency, 0.80)
    ),
    aes(label = Motif),
    size = 2.5,
    vjust = -1.5,
    check_overlap = TRUE
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # Reverse axis so hydrophobic motifs appear on the left
  scale_x_reverse(
    name = "Average Tripeptide Hydrophobicity (Hydrophobic ← → Hydrophilic)"
  ) +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_brewer(
    palette = "RdBu",
    name = "Net Charge"
  ) +
  
  labs(
    title = "Unordered Tripeptide Co-occurrence",
    subtitle = "Example: EPE, PEE and EEP are grouped together",
    x = NULL,
    y = "Frequency of Occurrence in Sliding Windows"
  ) +
  
  theme_minimal() +
  theme(legend.position = "right")

print(p_scatter_hydro)

# ------------------------------------------------------------
# PART 5: Save Output Figures
# ------------------------------------------------------------

# Vector format (ideal for publications)
ggsave(
  "results/tripeptide_enrichment.svg",
  plot = p_scatter_hydro,
  width = 8,
  height = 6
)

# Raster format (web/documents)
ggsave(
  "results/tripeptide_enrichment.png",
  plot = p_scatter_hydro,
  width = 8,
  height = 6,
  dpi = 300
)