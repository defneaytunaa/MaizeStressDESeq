# 0) Install & load necessary R/Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("clusterProfiler", "enrichplot", "org.Osativa.eg.db", "DOSE", "ggplot2")) 
# Use the organism-specific annotation package for maize; here: org.Osativa.eg.db (if available)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

#GO-Enrichment for reproductive system 
res <- read.csv("DESeq2_results_reproductive.csv", row.names=1)

# Define DEGs
deg_rep <- subset(res, padj < 0.05)

# Optional stronger filter:
deg_rep_strict <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

genes_rep <- rownames(deg_rep)
genes_rep_clean <- gsub("gene:", "", genes_rep)

install.packages("gprofiler2")
library(gprofiler2)

gost_rep <- gost(
  query = genes_rep_clean,
  organism = "zmays",
  correction_method = "fdr",
  sources = c("GO:BP","GO:MF","GO:CC")
)
gostplot(gost_rep, capped = TRUE)

gostres <- gost_rep


#Saving the results 
# Flatten the list columns in the gprofiler2 result table
gost_flat <- as.data.frame(
  lapply(gost_result, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) paste(x, collapse = ","))
    } else {
      col
    }
  })
)

write.csv(
  gost_flat,
  file = "GO_enrichment_reproductive_flat.csv",
  row.names = FALSE
)

# ==== same for vegetative system ==== 
res_veg <- read.csv("DESeq2_results_vegetative.csv",
                    header = TRUE,
                    row.names = 1)
head(res_veg)

deg_veg <- subset(res_veg, padj < 0.05)
genes_veg_raw <- rownames(deg_veg)
genes_veg <- gsub("gene:", "", genes_veg_raw)

library(gprofiler2)

gost_veg <- gost(
  query     = genes_veg,
  organism  = "zmays",
  correction_method = "fdr",
  sources   = c("GO:BP", "GO:MF", "GO:CC")
)
gostplot(gost_veg, capped = TRUE)

#save results
save_gost <- function(gost_obj, filename_prefix) {
  res <- gost_obj$result
  
  # 1) Flatten list columns so write.csv will work
  res_flat <- as.data.frame(
    lapply(res, function(col) {
      if (is.list(col)) {
        sapply(col, function(x) paste(x, collapse = ","))
      } else {
        col
      }
    })
  )
  
  # 2) Save full table (all columns)
  write.csv(
    res_flat,
    file = paste0(filename_prefix, "_full.csv"),
    row.names = FALSE
  )
  
  # 3) Save a simple, human-readable table
  #    Only use columns that actually exist
  wanted_cols <- c(
    "term_id",
    "term_name",
    "source",
    "p_value",           # adjusted p-value in gprofiler2
    "intersection_size",
    "query_size",
    "term_size",
    "precision",
    "recall"
  )
  
  simple_cols <- intersect(wanted_cols, colnames(res))
  simple <- res[, simple_cols]
  
  write.csv(
    simple,
    file = paste0(filename_prefix, "_simple.csv"),
    row.names = FALSE
  )
}

save_gost(gost_rep, "GO_enrichment_reproductive")
save_gost(gost_veg, "GO_enrichment_vegetative")