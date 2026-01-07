#========== GO Enrichment ==================

# strict DEGs already extracted earlier:
genes_40070 <- genes_40070_strict
genes_769   <- genes_769_strict

# Overlap (strict DEGs shared by both studies)
overlap_genes <- intersect(genes_40070, genes_769)

length(overlap_genes)
head(overlap_genes)

if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}
library(gprofiler2)

go_overlap <- gost(
  query = overlap_genes,
  organism = "zmays",
  correction_method = "fdr",
  sources = c("GO:BP", "GO:MF", "GO:CC")
)

gostplot(go_overlap, capped = TRUE)

# Full table
go_full <- go_overlap$result
write.csv(go_full, "GO_overlap_full.csv", row.names = FALSE)

# Simplified table
wanted_cols <- c(
  "term_id","term_name","source",
  "p_value","intersection_size",
  "term_size","query_size"
)

go_simple <- go_full[, intersect(wanted_cols, colnames(go_full))]
write.csv(go_simple, "GO_overlap_simple.csv", row.names = FALSE)

# ==== FLATTEN go_overlap$result BEFORE SAVING ====

go_full <- go_overlap$result

# Flatten list columns (gprofiler2 sometimes puts list-type columns)
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

# Save full flattened table
write.csv(go_full_flat, "GO_overlap_full.csv", row.names = FALSE)

# ==== Simple table ====

wanted_cols <- c(
  "term_id","term_name","source",
  "p_value","intersection_size",
  "term_size","query_size"
)

simple_cols <- intersect(wanted_cols, colnames(go_full_flat))

go_simple <- go_full_flat[, simple_cols, drop = FALSE]

write.csv(go_simple, "GO_overlap_simple.csv", row.names = FALSE)



## ==== KEGG enrichment for overlapping DEGs ====

# If not loaded yet:
if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}
library(gprofiler2)

kegg_overlap <- gost(
  query             = overlap_genes,
  organism          = "zmays",
  correction_method = "fdr",
  sources           = "KEGG"
)

# Quick visual check (gProfiler style)
gostplot(kegg_overlap, capped = TRUE)

## ==== Save KEGG overlap results ====

kegg_full <- kegg_overlap$result

# Flatten list columns so write.csv works
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

# Full table
write.csv(kegg_full_flat, "KEGG_overlap_full.csv", row.names = FALSE)

# Simplified table
wanted_cols <- c(
  "term_id", "term_name", "source",
  "p_value", "intersection_size",
  "term_size", "query_size"
)

simple_cols <- intersect(wanted_cols, colnames(kegg_full_flat))
kegg_simple <- kegg_full_flat[, simple_cols, drop = FALSE]

write.csv(kegg_simple, "KEGG_overlap_simple.csv", row.names = FALSE)

## ==== ggplot barplot of top KEGG overlap pathways ====

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Order by p-value, take top 10
kegg_simple$minus_log10_p <- -log10(kegg_simple$p_value)
kegg_simple <- kegg_simple[order(kegg_simple$p_value), ]
top_kegg <- head(kegg_simple, 10)

ggplot(top_kegg,
       aes(x = reorder(term_name, minus_log10_p),
           y = minus_log10_p)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top KEGG pathways – Overlapping drought DEGs (GSE40070 ∩ GSE76939)",
    x     = "KEGG pathway",
    y     = "-log10(adj p-value)"
  ) +
  theme_minimal(base_size = 12)

## ==== Bubble plot for KEGG overlap ====

top_kegg$GeneRatio <- top_kegg$intersection_size / top_kegg$term_size

ggplot(top_kegg,
       aes(x = GeneRatio,
           y = reorder(term_name, GeneRatio),
           size  = intersection_size,
           color = minus_log10_p)) +
  geom_point() +
  scale_color_gradient(low = "steelblue", high = "firebrick") +
  labs(
    title = "KEGG enrichment – Overlapping drought DEGs",
    x     = "Gene ratio (overlap DEGs / pathway genes)",
    y     = "KEGG pathway",
    size  = "DEG count",
    color = "-log10(adj p-value)"
  ) +
  theme_minimal(base_size = 12)
