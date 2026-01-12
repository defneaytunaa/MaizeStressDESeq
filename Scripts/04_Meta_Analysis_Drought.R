## ================================
## Helper: read DESeq2 and extract DEGs
## ================================

get_deg_sets <- function(file,
                         padj_cutoff   = 0.05,
                         log2fc_cutoff = 1) {
  stopifnot(file.exists(file))
  message("Reading: ", file)
  
  res <- read.csv(file,
                  header = TRUE,
                  row.names = 1,
                  check.names = FALSE)
  
  # keep only rows with valid padj
  res <- res[!is.na(res$padj), ]
  
  # up- and downregulated DEGs
  up_df   <- subset(res, padj < padj_cutoff & log2FoldChange >=  log2fc_cutoff)
  down_df <- subset(res, padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff)
  
  # function to extract gene IDs from a DESeq result data frame
  get_ids <- function(df) {
    if (nrow(df) == 0) return(character(0))
    
    # 1) prefer rownames
    if (!is.null(rownames(df)) && !all(rownames(df) == "")) {
      ids <- rownames(df)
    } else {
      # 2) fallback: first column that looks like an ID
      id_col <- grep("gene|id|Zm", colnames(df),
                     ignore.case = TRUE,
                     value = TRUE)[1]
      if (is.na(id_col)) {
        stop("Could not find a gene ID column in file: ", file)
      }
      ids <- as.character(df[[id_col]])
    }
    
    ids <- gsub("^gene:", "", ids)  # strip "gene:" prefix if present
    unique(ids)
  }
  
  up_ids   <- get_ids(up_df)
  down_ids <- get_ids(down_df)
  
  message("  Upregulated DEGs:   ", length(up_ids))
  message("  Downregulated DEGs: ", length(down_ids))
  
  list(
    up   = up_ids,
    down = down_ids
  )
}

## ================================
## Load DEGs for the two B73 drought datasets
## ================================

setwd("/Volumes/ThesisSSD/Results/GSE40070")  

deg_veg    <- get_deg_sets("DESeq2_results_vegetative.csv")

setwd("/Volumes/ThesisSSD/Results/GSE76939")  
res769 <- read.csv("DESeq2_Drought_vs_Control_FULL.csv", header = TRUE, check.names = FALSE)

deg_gse769 <- subset(res769, !is.na(padj) & padj < 0.05)

genes_gse769 <- gsub("gene:", "", deg_gse769$gene_id)

get_deg_sets <- function(file, padj_cutoff = 0.05, log2fc_cutoff = 1) {
  
  # Load without forcing rownames
  df <- read.csv(file, header = TRUE, check.names = FALSE)
  
  if (!"gene_id" %in% colnames(df)) stop("No gene_id column found in file!")
  
  # Basic DEG filter
  df_sig <- subset(df, !is.na(padj) & padj < padj_cutoff)
  
  # Optional strong filter
  df_strict <- subset(df_sig, abs(log2FoldChange) >= log2fc_cutoff)
  
  # Extract genes
  genes_all    <- gsub("gene:", "", df_sig$gene_id)
  genes_strict <- gsub("gene:", "", df_strict$gene_id)
  
  return(list(
    all_DEG_genes    = unique(genes_all),
    strict_DEG_genes = unique(genes_strict)
  ))
}

# All DEGs, ignoring direction
deg_veg_all    <- union(deg_veg$up,    deg_veg$down)
deg_gse769_all <- union(deg_gse769$up, deg_gse769$down)

## ===================================================
## 1) Compute gene-level overlaps (ALL DEGs)
## ===================================================

overlap_all <- intersect(deg_veg_all, deg_gse769_all)
unique_veg_all  <- setdiff(deg_veg_all,  deg_gse769_all)
unique_769_all  <- setdiff(deg_gse769_all, deg_veg_all)

cat("Overlap (all DEGs): ", length(overlap_all), "\n")
cat("Unique to GSE40070 (vegetative): ", length(unique_veg_all), "\n")
cat("Unique to GSE76939 (B73 drought): ", length(unique_769_all), "\n")

## ===================================================
## 2) Compute strict overlaps (padj<0.05 & |log2FC|≥1)
## ===================================================

overlap_strict <- intersect(deg_veg$strict_DEG_genes, deg_gse769$strict_DEG_genes)
unique_veg_strict <- setdiff(deg_veg$strict_DEG_genes, deg_gse769$strict_DEG_genes)
unique_769_strict <- setdiff(deg_gse769$strict_DEG_genes, deg_veg$strict_DEG_genes)

cat("Overlap (strict DEGs): ", length(overlap_strict), "\n")
cat("Unique strict GSE40070: ", length(unique_veg_strict), "\n")
cat("Unique strict GSE76939: ", length(unique_769_strict), "\n")

## ===================================================
## 3) Save overlap + unique lists for downstream GO/KEGG
## ===================================================

write.csv(overlap_all,       "Overlap_DEGs_ALL.csv",          row.names = FALSE)
write.csv(overlap_strict,    "Overlap_DEGs_STRICT.csv",       row.names = FALSE)
write.csv(unique_veg_all,    "Unique_GSE40070_ALL.csv",       row.names = FALSE)
write.csv(unique_769_all,    "Unique_GSE76939_ALL.csv",       row.names = FALSE)
write.csv(unique_veg_strict, "Unique_GSE40070_STRICT.csv",    row.names = FALSE)
write.csv(unique_769_strict, "Unique_GSE76939_STRICT.csv",    row.names = FALSE)


## ============================================
## Strict DEGs (padj < 0.05 & |log2FC| >= 1)
## for GSE40070 vegetative & GSE76939 B73 drought
## ============================================

library(VennDiagram)
library(grid)

## --- 1) Read GSE40070 vegetative DESeq2 results ---
setwd("/Volumes/ThesisSSD/Results/GSE40070")

res40070 <- read.csv(
  "DESeq2_results_vegetative.csv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)

res40070_strict <- subset(
  res40070,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
)

genes_40070_strict <- rownames(res40070_strict)
genes_40070_strict <- gsub("^gene:", "", genes_40070_strict)


## --- 2) Read GSE76939 drought vs control DESeq2 results ---
setwd("/Volumes/ThesisSSD/Results/GSE76939")

res769 <- read.csv(
  "DESeq2_Drought_vs_Control_FULL.csv",
  header      = TRUE,
  check.names = FALSE
)

res769_strict <- subset(
  res769,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
)

genes_769_strict <- gsub("^gene:", "", res769_strict$gene_id)


## --- 3) Venn diagram for STRICT DEGs ---
venn_strict <- venn.diagram(
  x = list(
    GSE40070 = genes_40070_strict,
    GSE76939 = genes_769_strict
  ),
  filename = NULL,
  fill     = c("skyblue", "firebrick"),
  alpha    = 0.5,
  cex      = 2,
  cat.cex  = 2,
  main     = "Overlap of STRICT DEGs: GSE40070 vs GSE76939"
)

grid.newpage()
grid.draw(venn_strict)


## =========================
## 1. GSE40070 – vegetative B73 drought vs control
## =========================

setwd("/Volumes/ThesisSSD/Results/GSE40070")

res_40070 <- read.csv(
  "DESeq2_results_vegetative.csv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)

# strict DEGs: padj < 0.05 & |log2FC| >= 1
deg_40070 <- subset(
  res_40070,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
)

genes_40070 <- gsub("^gene:", "", rownames(deg_40070))
genes_40070 <- unique(genes_40070)

cat("GSE40070 strict drought DEGs:", length(genes_40070), "\n")


## =========================
## 2. GSE76939 – B73 drought vs control
## =========================

setwd("/Volumes/ThesisSSD/Results/GSE76939")

res_769 <- read.csv(
  "DESeq2_Drought_vs_Control_FULL.csv",
  header      = TRUE,
  check.names = FALSE
)

# strict DEGs
deg_769 <- subset(
  res_769,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
)

# gene IDs are in the column 'gene_id'
genes_769 <- gsub("^gene:", "", deg_769$gene_id)
genes_769 <- unique(genes_769)

cat("GSE76939 strict drought DEGs:", length(genes_769), "\n")


## =========================
## 3. Compute overlap = drought core DEGs
## =========================

drought_core <- intersect(genes_40070, genes_769)
length(drought_core)
head(drought_core)

# build a small data frame for saving
drought_core_df <- data.frame(
  gene_id = sort(drought_core),
  stringsAsFactors = FALSE
)


## =========================
## 4. Save to Cross_Stress_Comparison folder
## =========================

setwd("/Volumes/ThesisSSD/Results/Study Comparison/Cross_Stress_Comparison")

write.csv(
  drought_core_df,
  file      = "Drought_core_DEGs.csv",
  row.names = FALSE
)

cat("Saved", nrow(drought_core_df),
    "drought core DEGs to Drought_core_DEGs.csv\n")
