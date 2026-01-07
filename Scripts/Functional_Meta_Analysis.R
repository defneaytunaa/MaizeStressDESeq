library(dplyr)

# helper: read a CSV that has one column with the maize IDs
read_gene_list <- function(file) {
  df <- read.csv(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Try to find a column that looks like IDs
  id_col <- grep("gene|id|Zm0000", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col)) {
    # if nothing matches, just take the first column
    id_col <- colnames(df)[1]
  }
  
  genes <- df[[id_col]]
  genes <- as.character(genes)
  genes <- gsub("^gene:", "", genes)   # remove "gene:" prefix if present
  genes <- trimws(genes)
  genes <- genes[genes != "" & !is.na(genes)]
  unique(genes)
}

genes_heat    <- read_gene_list("Heat_core_DEGs.csv")
genes_cold    <- read_gene_list("Cold_overlap_strict_genes.csv")
genes_drought <- read_gene_list("Drought_core_DEGs.csv")

length(genes_heat); length(genes_cold); length(genes_drought)
head(genes_heat)


if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}
library(gprofiler2)

run_go_bp <- function(genes, prefix) {
  if (length(genes) < 5) {
    message("Not enough genes for ", prefix, " – skipping.")
    return(NULL)
  }
  
  go_res <- gost(
    query             = genes,
    organism          = "zmays",
    correction_method = "fdr",
    sources           = c("GO:BP")   # only Biological Process
  )
  
  if (is.null(go_res) || is.null(go_res$result)) {
    message("No significant GO:BP enrichment for ", prefix)
    return(NULL)
  }
  
  res <- go_res$result
  
  # Flatten list columns so CSV works
  res_flat <- as.data.frame(
    lapply(res, function(col) {
      if (is.list(col)) sapply(col, function(x) paste(x, collapse = ",")) else col
    }),
    stringsAsFactors = FALSE
  )
  
  # Save full table
  write.csv(res_flat,
            file = paste0(prefix, "_GO_BP_full.csv"),
            row.names = FALSE)
  
  # Save a simpler table with just key columns
  keep_cols <- c("term_id", "term_name", "p_value",
                 "intersection_size", "query_size", "term_size")
  keep_cols <- intersect(keep_cols, colnames(res_flat))
  
  simple <- res_flat[, keep_cols, drop = FALSE]
  write.csv(simple,
            file = paste0(prefix, "_GO_BP_simple.csv"),
            row.names = FALSE)
  
  message("Saved GO:BP results for ", prefix,
          " (", nrow(simple), " terms).")
  invisible(simple)
}

go_heat_simple    <- run_go_bp(genes_heat,    "Heat_core")
go_cold_simple    <- run_go_bp(genes_cold,    "Cold_overlap")
go_drought_simple <- run_go_bp(genes_drought, "Drought_core")

#=== as go_heat_simple did not have any significant results, I will continue with cold and drought" 
library(gprofiler2)

run_go_bp <- function(genes, prefix, significant = TRUE) {
  if (length(genes) < 5) {
    message("Not enough genes for ", prefix, " – skipping.")
    return(NULL)
  }
  
  message("Running GO:BP for ", prefix,
          " (significant = ", significant, ")")
  
  go_res <- gost(
    query             = genes,
    organism          = "zmays",
    correction_method = "fdr",
    sources           = c("GO:BP"),
    significant       = significant   # <--- key difference
  )
  
  # Check if we actually got results
  if (is.null(go_res) || is.null(go_res$result) || nrow(go_res$result) == 0) {
    message("No GO:BP enrichment results for ", prefix,
            " under these settings.")
    return(NULL)
  }
  
  res <- go_res$result
  
  # Flatten list columns for saving
  res_flat <- as.data.frame(
    lapply(res, function(col) {
      if (is.list(col)) sapply(col, function(x) paste(x, collapse = ",")) else col
    }),
    stringsAsFactors = FALSE
  )
  
  # Save full table
  full_file <- paste0(prefix, "_GO_BP_full.csv")
  write.csv(res_flat, full_file, row.names = FALSE)
  
  # Save simpler table
  keep_cols <- c("term_id", "term_name", "p_value",
                 "intersection_size", "query_size", "term_size")
  keep_cols <- intersect(keep_cols, colnames(res_flat))
  simple <- res_flat[, keep_cols, drop = FALSE]
  
  simple_file <- paste0(prefix, "_GO_BP_simple.csv")
  write.csv(simple, simple_file, row.names = FALSE)
  
  message("Saved ", nrow(simple), " GO:BP terms for ", prefix,
          " to: ", simple_file)
  
  invisible(simple)
}

# Heat: exploratory (no FDR filter)
go_heat_simple    <- run_go_bp(genes_heat,    "Heat_core",    significant = FALSE)

# Cold & Drought: proper FDR-significant enrichment
go_cold_simple    <- run_go_bp(genes_cold,    "Cold_overlap", significant = TRUE)
go_drought_simple <- run_go_bp(genes_drought, "Drought_core", significant = TRUE)

#== Combining Results ==
library(dplyr)
library(tidyr)

go_list <- list(
  Heat    = go_heat_simple,
  Cold    = go_cold_simple,
  Drought = go_drought_simple
)

# Bind only non-NULL tables and tag them with stress
go_all <- bind_rows(
  lapply(names(go_list), function(stress_name) {
    df <- go_list[[stress_name]]
    if (is.null(df)) return(NULL)
    df$stress <- stress_name
    df
  })
)

# Build meta-matrix: -log10(p_value) per stress
go_wide <- go_all %>%
  mutate(neglog10p = -log10(p_value)) %>%
  select(term_id, term_name, stress, neglog10p) %>%
  pivot_wider(
    id_cols = c(term_id, term_name),
    names_from  = stress,
    values_from = neglog10p,
    values_fill = 0
  )

write.csv(go_wide, "GO_BP_meta_matrix_heat_cold_drought.csv", row.names = FALSE)
head(go_wide)

## =============================
## Run GO:BP enrichment (functional meta-analysis)
## =============================

# HEAT — exploratory mode (no FDR significance filtering)
go_heat_simple <- run_go_bp(
  genes_heat,
  prefix      = "Heat_core",
  significant = FALSE
)

# COLD — real FDR-significant enrichment
go_cold_simple <- run_go_bp(
  genes_cold,
  prefix      = "Cold_overlap",
  significant = TRUE
)

# DROUGHT — real FDR-significant enrichment
go_drought_simple <- run_go_bp(
  genes_drought,
  prefix      = "Drought_core",
  significant = TRUE
)


## =============================
## 0) Libraries
## =============================
library(dplyr)
library(tidyr)

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

## =============================
## 1) Set working directory
## =============================
setwd("/Volumes/ThesisSSD/Results/Cross_Stress_Comparison")


## =============================
## 2) Load GO:BP simple results
## =============================
go_cold    <- read.csv("Cold_overlap_GO_BP_simple.csv", stringsAsFactors = FALSE)
go_drought <- read.csv("Drought_core_GO_BP_simple.csv", stringsAsFactors = FALSE)

go_cold$stress    <- "Cold"
go_drought$stress <- "Drought"


## =============================
## 3) Combine into one table
## =============================
go_all <- bind_rows(go_cold, go_drought)

## Convert to –log10(p)
go_all <- go_all %>%
  mutate(neglog10p = -log10(p_value))

## =============================
## 4) Build wide matrix
## =============================
go_wide <- go_all %>%
  select(term_id, term_name, stress, neglog10p) %>%
  pivot_wider(
    id_cols    = c(term_id, term_name),
    names_from = stress,
    values_from = neglog10p,
    values_fill = 0
  )

# Save for reference
write.csv(go_wide, "GO_BP_meta_matrix_Heat_Cold_Drought.csv",
          row.names = FALSE)


## =============================
## 5) Extract numeric matrix
## =============================
mat <- as.matrix(go_wide[, c("Cold", "Drought")])
rownames(mat) <- go_wide$term_name


## =============================
## 6) Select top N most enriched terms
## =============================
topN <- 50   # change if you want 20 / 100 etc.

# rank rows by maximum enrichment value
sel_idx <- order(apply(mat, 1, max), decreasing = TRUE)[1:topN]

mat_sel <- mat[sel_idx, ]  # <-- THIS is the object needed for the heatmap


## =============================
## 7) Plot heatmap
## =============================
pheatmap(
  mat_sel,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "GO:BP Functional Signature – Cold vs Drought",
  fontsize_row = 6,
  fontsize_col = 10,
  angle_col = 0,
  border_color = NA
)