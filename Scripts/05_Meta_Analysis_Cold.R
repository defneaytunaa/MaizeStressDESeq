## ============================================
## Cold stress – strict DEGs & overlap (B73)
## ============================================

library(VennDiagram)
library(grid)

## --- 1) GSE225916: B73 cold vs control ---

setwd("/Volumes/ThesisSSD/Results/GSE225916")  

res_cold_225916 <- read.csv(
  "DESeq2_results_Cold_vs_Control_FULL.csv",   # adjust name if different
  header      = TRUE,
  check.names = FALSE
)

res_cold_225916_strict <- subset(
  res_cold_225916,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
)

# gene IDs (column 'gene_id' or whatever you used)
genes_225916_cold <- gsub("^gene:", "", res_cold_225916_strict$gene_id)


## --- 2) GSE76939: B73 cold vs control ---

setwd("/Volumes/ThesisSSD/Results/GSE76939")    # adjust path if needed

res_cold_76939 <- read.csv(
  "DESeq2_results_sourceName_Cold_vs_Control_full.csv",  # adjust name if needed
  header      = TRUE,
  check.names = FALSE
)

res_cold_76939_strict <- subset(
  res_cold_76939,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
)

genes_76939_cold <- gsub("^gene:", "", res_cold_76939_strict$gene_id)


## --- 3) Overlap of strict cold DEGs ---

cold_overlap_genes <- intersect(genes_225916_cold, genes_76939_cold)
cold_only_225916   <- setdiff(genes_225916_cold, genes_76939_cold)
cold_only_76939    <- setdiff(genes_76939_cold, genes_225916_cold)

cat("Strict cold DEGs GSE225916: ", length(genes_225916_cold), "\n")
cat("Strict cold DEGs GSE76939:  ", length(genes_76939_cold), "\n")
cat("Overlap (strict) cold DEGs: ", length(cold_overlap_genes), "\n")


## --- 4) Save the overlap & unique lists ---

write.csv(data.frame(gene_id = cold_overlap_genes),
          "Cold_overlap_strict_genes.csv",
          row.names = FALSE)

write.csv(data.frame(gene_id = cold_only_225916),
          "Cold_unique_GSE225916_strict.csv",
          row.names = FALSE)

write.csv(data.frame(gene_id = cold_only_76939),
          "Cold_unique_GSE76939_strict.csv",
          row.names = FALSE)


## --- 5) Venn diagram – strict cold DEGs ---

venn_cold <- venn.diagram(
  x = list(
    GSE225916 = genes_225916_cold,
    GSE76939  = genes_76939_cold
  ),
  filename = NULL,
  fill     = c("skyblue", "tomato"),
  alpha    = 0.5,
  cex      = 2,
  cat.cex  = 2,
  main     = "Overlap of STRICT DEGs – Cold vs Control (B73)"
)

grid.newpage()
grid.draw(venn_cold)



## ============================================
## GO enrichment – overlapping cold DEGs
## ============================================

if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}
library(gprofiler2)

# Make sure the overlap is defined
length(cold_overlap_genes)
head(cold_overlap_genes)

go_cold_overlap <- gost(
  query             = cold_overlap_genes,
  organism          = "zmays",
  correction_method = "fdr",
  sources           = c("GO:BP", "GO:MF", "GO:CC")
)

# Quick visual
gostplot(go_cold_overlap, capped = TRUE)


## ---- Save GO results (flatten lists) ----

go_full <- go_cold_overlap$result

go_full_flat <- as.data.frame(
  lapply(go_full, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) paste(x, collapse = ","))
    } else {
      col
    }
  }),
  stringsAsFactors = FALSE
)

write.csv(go_full_flat, "GO_Cold_overlap_full.csv", row.names = FALSE)

wanted_cols <- c(
  "term_id","term_name","source",
  "p_value","intersection_size",
  "term_size","query_size"
)

simple_cols <- intersect(wanted_cols, colnames(go_full_flat))
go_cold_simple <- go_full_flat[, simple_cols, drop = FALSE]

write.csv(go_cold_simple, "GO_Cold_overlap_simple.csv", row.names = FALSE)




## ============================================
## KEGG enrichment – overlapping cold DEGs
## ============================================

kegg_cold_overlap <- gost(
  query             = cold_overlap_genes,
  organism          = "zmays",
  correction_method = "fdr",
  sources           = "KEGG"
)

gostplot(kegg_cold_overlap, capped = TRUE)


## ---- Save KEGG results (flatten lists) ----

kegg_full <- kegg_cold_overlap$result

kegg_full_flat <- as.data.frame(
  lapply(kegg_full, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) paste(x, collapse = ","))
    } else {
      col
    }
  }),
  stringsAsFactors = FALSE
)

write.csv(kegg_full_flat, "KEGG_Cold_overlap_full.csv", row.names = FALSE)

wanted_cols <- c(
  "term_id","term_name","source",
  "p_value","intersection_size",
  "term_size","query_size"
)

simple_cols <- intersect(wanted_cols, colnames(kegg_full_flat))
kegg_cold_simple <- kegg_full_flat[, simple_cols, drop = FALSE]

write.csv(kegg_cold_simple, "KEGG_Cold_overlap_simple.csv", row.names = FALSE)



## ============================================
## Dotplot-style KEGG plot – Cold overlap
## ============================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Make sure kegg_cold_simple exists
# If needed, re-read:
# kegg_cold_simple <- read.csv("KEGG_Cold_overlap_simple.csv")

# Add -log10(p) and GeneRatio
kegg_cold_simple$minus_log10_p <- -log10(kegg_cold_simple$p_value)
kegg_cold_simple$GeneRatio <- kegg_cold_simple$intersection_size / kegg_cold_simple$term_size

# Order by significance (lowest p -> top)
kegg_cold_simple <- kegg_cold_simple[order(kegg_cold_simple$p_value), ]

# Pick top N pathways to avoid clutter
top_n <- 20
top_kegg_cold <- head(kegg_cold_simple, top_n)

ggplot(top_kegg_cold,
       aes(x = GeneRatio,
           y = reorder(term_name, GeneRatio),
           size  = intersection_size,
           color = minus_log10_p)) +
  geom_point() +
  scale_color_gradient(low = "steelblue", high = "firebrick") +
  labs(
    title = "KEGG enrichment – Overlapping cold DEGs (B73)",
    x     = "Gene ratio (overlap DEGs / pathway genes)",
    y     = "KEGG pathway",
    size  = "DEG count (overlap)",
    color = "-log10(adj p-value)"
  ) +
  theme_minimal(base_size = 12)