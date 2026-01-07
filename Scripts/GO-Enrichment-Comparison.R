
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("jrwalsh/MaizeGO")
library(MaizeGO)
data("MaizeGO.B73.v4") 

## --- Working directory with your three lists ---
setwd("/Volumes/ThesisSSD/Results/Study_Comparison/Cross_Stress_Comparison")

# Check what's there
list.files()

# 1) Read gene lists (assume first column has gene IDs)
heat_df    <- read.csv("Heat_core_DEGs.csv",              header = TRUE, stringsAsFactors = FALSE)
cold_df    <- read.csv("Cold_overlap_strict_genes.csv",   header = TRUE, stringsAsFactors = FALSE)
drought_df <- read.csv("Drought_core_DEGs.csv",           header = TRUE, stringsAsFactors = FALSE)

genes_heat    <- unique(as.character(heat_df[,1]))
genes_cold    <- unique(as.character(cold_df[,1]))
genes_drought <- unique(as.character(drought_df[,1]))

length(genes_heat); length(genes_cold); length(genes_drought)

## --- MaizeGO annotation ---
library(MaizeGO)
data("MaizeGO.B73.v4")
head(MaizeGO.B73.v4)
colnames(MaizeGO.B73.v4)

library(dplyr)

# Quick sanity check
colnames(MaizeGO.B73.v4)
unique(MaizeGO.B73.v4$type)

# Keep only GO:BP rows
go_bp <- MaizeGO.B73.v4 %>%
  dplyr::filter(type == "BP")

library(dplyr)
library(clusterProfiler)

# Ensure gene IDs are comparable (remove prefixes if needed)
clean_id <- function(x) gsub("^gene:", "", x)

genes_heat     <- clean_id(genes_heat)
genes_cold     <- clean_id(genes_cold)
genes_drought  <- clean_id(genes_drought)

# ---- Map genes --> GO terms ----
map_to_go <- function(gene_list, go_table) {
  go_table %>%
    dplyr::filter(geneID %in% gene_list) %>%
    dplyr::select(geneID, goTerm) %>%
    distinct()
}

go_heat     <- map_to_go(genes_heat, go_bp)
go_cold     <- map_to_go(genes_cold, go_bp)
go_drought  <- map_to_go(genes_drought, go_bp)

# Quick sanity check:
nrow(go_heat); nrow(go_cold); nrow(go_drought)

go_terms_heat    <- unique(go_heat$goTerm)
go_terms_cold    <- unique(go_cold$goTerm)
go_terms_drought <- unique(go_drought$goTerm)

library(GOSemSim)

# Build a custom GO annotation object from your maize file
maize_annot <- go_bp %>%
  dplyr::select(geneID, goTerm, type)

# Build the semData object
maize_sem <- godata(
  ont = "BP",
  keytype = "GENEID",
  organism = maize_annot,
  computeIC = TRUE
)

go_heat_terms     <- go_heat_simple$term_id
go_cold_terms     <- go_cold_simple$term_id
go_drought_terms  <- go_drought_simple$term_id

all_terms <- unique(c(go_heat_terms, go_cold_terms, go_drought_terms))

mat <- data.frame(
  Heat    = as.integer(all_terms %in% go_heat_terms),
  Cold    = as.integer(all_terms %in% go_cold_terms),
  Drought = as.integer(all_terms %in% go_drought_terms),
  row.names = all_terms
)

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

sim_heat_cold    <- jaccard(go_heat_terms, go_cold_terms)
sim_heat_drought <- jaccard(go_heat_terms, go_drought_terms)
sim_cold_drought <- jaccard(go_cold_terms, go_drought_terms)

sim_matrix <- matrix(
  c(1, sim_heat_cold, sim_heat_drought,
    sim_heat_cold, 1, sim_cold_drought,
    sim_heat_drought, sim_cold_drought, 1),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("Heat", "Cold", "Drought"),
                  c("Heat", "Cold", "Drought"))
)

library(pheatmap)

pheatmap(
  sim_matrix,
  main = "Functional Similarity (GO:BP Enrichment)\nHeat vs Cold vs Drought",
  display_numbers = TRUE,
  fontsize = 12,
  color = colorRampPalette(c("white","steelblue"))(100)
)
