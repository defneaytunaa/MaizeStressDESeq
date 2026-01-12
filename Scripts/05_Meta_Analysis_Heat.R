## ============================================
## Build HEAT core DEG list (CM1 ∩ CB25 ∩ F1)
## using existing DESeq2 result files
## ============================================

# 1) Set working directory where your DESeq2 results for GSE122866 live
setwd("/Volumes/ThesisSSD/Results/GSE122866")  # <-- adjust if needed


# 2) Helper: read one DESeq2 result file and return strict DEGs (gene IDs)
get_strict_degs <- function(file,
                            padj_cutoff   = 0.05,
                            log2fc_cutoff = 1) {
  stopifnot(file.exists(file))
  message("Reading: ", file)
  
  df <- read.csv(file, header = TRUE, check.names = FALSE)
  
  # Filter by padj and |log2FC|
  if (!all(c("padj","log2FoldChange") %in% colnames(df))) {
    stop("File ", file, " does not contain 'padj' and 'log2FoldChange' columns.")
  }
  df_strict <- subset(df, !is.na(padj) & padj < padj_cutoff &
                        abs(log2FoldChange) >= log2fc_cutoff)
  
  if (nrow(df_strict) == 0) {
    warning("No strict DEGs found in ", file)
    return(character(0))
  }
  
  # Try to get gene IDs:
  # 1) If a column called 'gene_id' exists, use that
  # 2) Otherwise fall back to rownames
  if ("gene_id" %in% colnames(df_strict)) {
    ids <- as.character(df_strict$gene_id)
  } else {
    # if row.names were saved as a column during write.csv,
    # you might have something like 'X' or an unnamed first column
    # but we FIRST try rownames:
    rn <- rownames(df_strict)
    if (!is.null(rn) && !all(rn == "")) {
      ids <- rn
    } else {
      # fallback: use first column
      ids <- as.character(df_strict[[1]])
      warning("Using first column as gene IDs for file: ", file)
    }
  }
  
  # Clean "gene:" prefix if present
  ids <- gsub("^gene:", "", ids)
  ids <- unique(ids)
  
  message("  Strict DEGs found: ", length(ids))
  return(ids)
}


# 3) Read strict DEGs for each genotype
#    Adjust these filenames to match your actual DESeq2 result files
deg_cm1  <- get_strict_degs("DESeq2_CM1_heat_vs_control_FULL.csv")
deg_cb25 <- get_strict_degs("DESeq2_CB25_heat_vs_control_FULL.csv")
deg_f1   <- get_strict_degs("DESeq2_F1_heat_vs_control_FULL.csv")

length(deg_cm1)
length(deg_cb25)
length(deg_f1)


# 4) Compute HEAT CORE DEGs = intersection of all three genotypes
heat_core <- Reduce(intersect, list(deg_cm1, deg_cb25, deg_f1))

cat("HEAT core DEGs (CM1 ∩ CB25 ∩ F1): ", length(heat_core), "\n")


# 5) (Optional) Also compute union & genotype–specific sets for later use
heat_union        <- unique(c(deg_cm1, deg_cb25, deg_f1))
heat_cm1_specific <- setdiff(deg_cm1,  union(deg_cb25, deg_f1))
heat_cb25_specific<- setdiff(deg_cb25, union(deg_cm1,  deg_f1))
heat_f1_specific  <- setdiff(deg_f1,   union(deg_cm1,  deg_cb25))


# 6) Save the core and supporting lists for meta-analysis
write.csv(data.frame(gene_id = heat_core),
          "Heat_core_DEGs.csv",
          row.names = FALSE)

write.csv(data.frame(gene_id = heat_union),
          "Heat_union_DEGs.csv",
          row.names = FALSE)

write.csv(data.frame(gene_id = heat_cm1_specific),
          "Heat_CM1_specific_DEGs.csv",
          row.names = FALSE)

write.csv(data.frame(gene_id = heat_cb25_specific),
          "Heat_CB25_specific_DEGs.csv",
          row.names = FALSE)

write.csv(data.frame(gene_id = heat_f1_specific),
          "Heat_F1_specific_DEGs.csv",
          row.names = FALSE)



## =========================================================
## 1) Load DEG lists (Heat, Cold, Drought)
## =========================================================

cold_df    <- read.csv("Cold_overlap_strict_genes.csv",    header = TRUE)
drought_df <- read.csv("Drought_core_DEGs.csv",             header = TRUE)
heat_df    <- read.csv("Heat_core_DEGs.csv",                header = TRUE)

cold_genes    <- unique(cold_df$gene_id)
drought_genes <- unique(drought_df$gene_id)
heat_genes    <- unique(heat_df$gene_id)

length(cold_genes)
length(drought_genes)
length(heat_genes)


## =========================================================
## 2) Venn Diagram (Cold vs Drought vs Heat)
## =========================================================

library(VennDiagram)
library(grid)

venn <- venn.diagram(
  x = list(
    Cold    = cold_genes,
    Drought = drought_genes,
    Heat    = heat_genes
  ),
  filename = NULL,
  fill     = c("skyblue", "gold", "tomato"),
  alpha    = 0.5,
  cex      = 2,
  cat.cex  = 2,
  main     = "Cross-Stress DEG Overlap (Cold / Drought / Heat)"
)

grid.newpage()
grid.draw(venn)


## =========================================================
## 3) Universal Abiotic Stress Core Genes
## (Genes shared in all 3 stress types)
## =========================================================

core_3way <- Reduce(intersect, list(cold_genes, drought_genes, heat_genes))
length(core_3way)

write.csv(
  data.frame(gene_id = core_3way),
  "StressCore_DEGs_Cold_Drought_Heat.csv",
  row.names = FALSE
)


## =========================================================
## 4) Pairwise overlaps
## =========================================================

cold_drought <- intersect(cold_genes, drought_genes)
cold_heat    <- intersect(cold_genes, heat_genes)
drought_heat <- intersect(drought_genes, heat_genes)

write.csv(data.frame(gene_id = cold_drought),
          "Overlap_Cold_Drought.csv", row.names = FALSE)

write.csv(data.frame(gene_id = cold_heat),
          "Overlap_Cold_Heat.csv", row.names = FALSE)

write.csv(data.frame(gene_id = drought_heat),
          "Overlap_Drought_Heat.csv", row.names = FALSE)


## =========================================================
## 5) Stress-Specific DEGs
## =========================================================

cold_specific    <- setdiff(cold_genes,    union(drought_genes, heat_genes))
drought_specific <- setdiff(drought_genes, union(cold_genes,    heat_genes))
heat_specific    <- setdiff(heat_genes,    union(cold_genes,    drought_genes))

write.csv(data.frame(gene_id = cold_specific),
          "Cold_specific_DEGs.csv", row.names = FALSE)

write.csv(data.frame(gene_id = drought_specific),
          "Drought_specific_DEGs.csv", row.names = FALSE)

write.csv(data.frame(gene_id = heat_specific),
          "Heat_specific_DEGs2.csv", row.names = FALSE)   # renamed to avoid confusion


## ================================================================
## Load the three DEG sets
## ================================================================
cold  <- read.csv("Cold_overlap_strict_genes.csv", header = TRUE, stringsAsFactors = FALSE)
drought <- read.csv("Drought_core_DEGs.csv", header = TRUE, stringsAsFactors = FALSE)
heat <- read.csv("Heat_core_DEGs.csv", header = TRUE, stringsAsFactors = FALSE)

## Each file likely has a column containing gene IDs.
## We detect that automatically:
get_ids <- function(df) {
  # find column containing gene IDs
  id_col <- grep("gene|id|Zm", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col)) stop("Could not find gene ID column!")
  
  ids <- df[[id_col]]
  ids <- gsub("^gene:", "", ids)   # clean gene: prefix if present
  unique(ids)
}

cold_genes    <- get_ids(cold)
drought_genes <- get_ids(drought)
heat_genes    <- get_ids(heat)

## ================================================================
## 1) Identify overlapping genes between Cold and Drought
## ================================================================
cold_drought_overlap <- intersect(cold_genes, drought_genes)
length(cold_drought_overlap)

## ================================================================
## 2) Exclude any gene that is ALSO a heat DEG (optional but recommended)
## ================================================================
cold_drought_specific <- setdiff(cold_drought_overlap, heat_genes)
length(cold_drought_specific)

## ================================================================
## 3) Save the final list
## ================================================================
write.csv(
  data.frame(gene_id = cold_drought_specific),
  "Cold_Drought_shared_DEGs.csv",
  row.names = FALSE
)