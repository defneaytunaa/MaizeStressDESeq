# ========= DESeq2 Script ========== 

# First convert the Data to CSV
# set working direcrory 
setwd("/Volumes/ThesisSSD/Results/GSE40070")
counts_raw <- read_csv("GSE40070_salmon.merged.gene_counts.csv")

library(readr)
library(dplyr)
library(DESeq2)

# Step 1 - Re-import with semicolon delimiter
counts_raw <- read_delim(
  "salmon.merged.gene_counts.tsv",
  delim = "\t",
  trim_ws = TRUE
)

# Step 2 - Inspect
print(colnames(counts_raw)[1:10])
print(dim(counts_raw))

# Step 3 — Build count_mat (drop gene_name, keep gene_id as rownames)
library(dplyr)
library(DESeq2)

# Keep IDs, drop annotation col(s), keep only the six sample columns
id_cols   <- c("gene_id","gene_name")
stopifnot(all(id_cols[1] %in% colnames(counts_raw)))  # gene_id must exist

gene_ids  <- counts_raw[["gene_id"]]
count_df  <- counts_raw %>% select(-any_of(id_cols))

# Safety: ensure all remaining columns are numeric integers (DESeq2 requires integers)
count_df[] <- lapply(count_df, function(x) {
  y <- suppressWarnings(as.numeric(x))
  if (any(is.na(y) & !is.na(x))) stop("Non-numeric entries found in count columns.")
  as.integer(round(y))
})

# Build matrix with gene IDs as rownames
count_mat <- as.matrix(count_df)
rownames(count_mat) <- make.unique(as.character(gene_ids))

# Quick checks
dim(count_mat)           # should be 43093 x 6
head(colnames(count_mat))
head(rownames(count_mat))
summary(as.vector(count_mat))

# === To sucessfully conduct the DESeq2 workflow, you need to define the conditions -> for that download the Metadata in the SRA Toolkit (convert to CSV first) ===
# Preprocessing of the metadata 
# --- Locate the metadata file robustly ---
cands <- c("SraRunTable.csv", "SraRunTable")
meta_path <- cands[file.exists(cands)]
if (length(meta_path) == 0) stop("Couldn't find SraRunTable(.csv) in the current folder.")
meta_path <- meta_path[1]
message("Reading: ", meta_path)

# --- Read metadata lines (Numbers exports often have a title line first) ---
meta <- read.table(
  meta_path,
  sep = ";",
  header=TRUE
)


# Lets split the data into the two different analysis
print(meta[,c("Run", "source_name")])
reproductive_samples <- c("SRR536834", "SRR536835", "SRR536836", "SRR536837")
vegetative_samples <- c("SRR536838", "SRR536839", "SRR536840", "SRR536841")

count_mat_reproductive <- count_mat[,reproductive_samples]
count_mat_vegetative <- count_mat[,vegetative_samples]

meta_reproductive <- meta[meta$Run %in% reproductive_samples,] 
meta_vegetative <- meta[meta$Run %in% vegetative_samples,]

# Lets first do it for the reproductive samples
# === now we can wire the conditions into the DESeq2 
# === Step 1: align metadata to your count matrix and build colData ==== 
# columns (samples) in your count matrix:
samples <- colnames(count_mat_reproductive)
samples
# Should be: "SRR536834" "SRR536835" "SRR536836" "SRR536839" "SRR536840" "SRR536841"

# Use the 'growth_condition' column as the grouping factor
stopifnot("growth_condition" %in% colnames(meta_reproductive))
condition <- factor(trimws(meta_reproductive$growth_condition))

# quick sanity check
table(condition, useNA="ifany")
levels(condition)

# Build colData
coldata <- data.frame(row.names = samples, condition = condition)
# Make sure reference level is "well watered"
coldata$condition <- relevel(coldata$condition, ref = "well watered")
coldata

# === Step B — build DESeqDataSet and run DESeq2 === 
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat_reproductive,
  colData   = coldata,
  design    = ~ condition
)

# optional prefilter (recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

# Contrast: drought stress vs well watered
res <- results(dds, contrast = c("condition","drought stress","well watered"))
summary(res)

# Peek at top DE genes (sorted by adjusted p-value)
head(res[order(res$padj), ], 10)

# ==== Save results: 
res_ordered <- res[order(res$padj), ]
write.csv(as.data.frame(res_ordered), "DESeq2_results_reproductive.csv")

sig <- subset(as.data.frame(res_ordered), padj < 0.05 & abs(log2FoldChange) >= 1)
write.csv(sig, "DESeq2_results_sig.csv")

nrow(sig)

# Visualisation with MA-Plot 
plotMA(res, ylim = c(-5, 5), main = "Drought stress vs Well watered") # -> this might be not concise enough 

 plotMA(
       res,
       ylim = c(-3, 3),             
       alpha = 0.1,                
       main = "Drought stress vs Well watered (padj < 0.1)"
   )
sig <- subset(res_df, padj < 0.05)

# for the shrinking:
# See all available coefficients
resultsNames(dds)
coef_name <- grep("condition.*drought.*vs.*well", resultsNames(dds), value = TRUE)
stopifnot(length(coef_name) == 1)

library(apeglm)
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

plotMA(res_shrunk, ylim = c(-6, 6), alpha = 0.05,
       main = "MA (apeglm-shrunk)")