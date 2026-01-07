# Install if missing
if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}

library(gprofiler2)

# Load overlapping genes
overlap <- read.csv("Cold_Drought_shared_DEGs.csv", stringsAsFactors = FALSE)

# Detect gene column automatically
id_col <- grep("gene|id|Zm", colnames(overlap), ignore.case = TRUE, value = TRUE)[1]
genes <- overlap[[id_col]]
genes <- gsub("^gene:", "", genes)   # clean prefix if present

# Run GO enrichment
go_cd <- gost(
  query            = genes,
  organism         = "zmays",
  correction_method = "fdr",
  sources          = c("GO:BP", "GO:MF", "GO:CC")
)

# Plot
gostplot(go_cd, capped = TRUE)

# Save results
go_full <- go_cd$result

# Flatten list columns
go_flat <- as.data.frame(
  lapply(go_full, function(col) {
    if (is.list(col)) sapply(col, paste, collapse = ",") else col
  }),
  stringsAsFactors = FALSE
)

write.csv(go_flat, "GO_Cold_Drought_overlap_full.csv", row.names = FALSE)

# Save simple version
wanted_cols <- c("term_id","term_name","source","p_value","intersection_size","term_size")
go_simple <- go_flat[, intersect(wanted_cols, colnames(go_flat))]
write.csv(go_simple, "GO_Cold_Drought_overlap_simple.csv", row.names = FALSE)


#===== KEGG Enruchment =======
library(gprofiler2)

# 1. Load overlapping genes
overlap <- read.csv("Cold_Drought_shared_DEGs.csv", stringsAsFactors = FALSE)

id_col <- grep("gene|id|Zm", colnames(overlap), ignore.case = TRUE, value = TRUE)[1]
genes  <- overlap[[id_col]]
genes  <- gsub("^gene:", "", genes)
genes  <- unique(genes)

length(genes)
head(genes)

kegg_cd <- gost(
  query             = genes,
  organism          = "zmays",
  correction_method = "fdr",
  sources           = "KEGG"
)

str(kegg_cd)