install.packages("gprofiler2")
library(gprofiler2)

genes_rep <- gsub("gene:", "", rownames(deg_rep))
genes_veg <- gsub("gene:", "", rownames(deg_veg))

#Run KEGG Enrichment 
kegg_rep <- gost(
  query             = genes_rep,
  organism          = "zmays",
  correction_method = "fdr",
  sources           = "KEGG"
)

kegg_veg <- gost(
  query             = genes_veg,
  organism          = "zmays",
  correction_method = "fdr",
  sources           = "KEGG"
)

#check results
str(kegg_rep$result)
str(kegg_veg$result)

# If non-null:
head(kegg_rep$result)
head(kegg_veg$result)

gostplot(kegg_rep, capped = TRUE)   # Reproductive
gostplot(kegg_veg, capped = TRUE)   # Vegetative

#save results 
save_kegg <- function(gost_obj, prefix) {
  res <- gost_obj$result
  if (is.null(res)) {
    message("No significant KEGG enrichment for ", prefix)
    return(NULL)
  }
  res_flat <- as.data.frame(lapply(res, function(col) {
    if (is.list(col)) sapply(col, paste, collapse=",") else col
  }))
  write.csv(res_flat, paste0(prefix, "_kegg_full.csv"), row.names = FALSE)
}
save_kegg(kegg_rep, "reproductive")
save_kegg(kegg_veg,   "vegetative")

#Barplot for better visuals 
library(ggplot2)

# Load KEGG results
kegg_rep <- read.csv("vegetative_kegg_full.csv")

# Select top pathways (e.g., top 10)
top_rep <- kegg_rep[order(kegg_rep$p_value), ][1:10, ]

# Plot
ggplot(top_rep, aes(x = reorder(term_name, -log10(p_value)), 
                    y = -log10(p_value))) +
  geom_col(fill = "tomato") +
  coord_flip() +
  labs(
    title = "Top KEGG Pathways – Vegetative System",
    x = "KEGG Pathway",
    y = "-log10(adj p-value)"
  ) +
  theme_minimal()




library(ggplot2)
veg_kegg <- read.csv("vegetative_kegg_full.csv")
veg_kegg <- veg_kegg[order(veg_kegg$p_value), ]
top_veg  <- head(veg_kegg, 10)

ggplot(top_veg,
       aes(x = reorder(term_name, -log10(p_value)),
           y = -log10(p_value))) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top KEGG pathways – Vegetative system",
    x     = "KEGG pathway",
    y     = "-log10(adj p-value)"
  ) +
  theme_minimal(base_size = 12)

top_veg$GeneRatio <- top_veg$intersection_size / top_veg$term_size

ggplot(top_veg,
       aes(x = GeneRatio,
           y = reorder(term_name, GeneRatio),
           size  = intersection_size,
           color = -log10(p_value))) +
  geom_point() +
  scale_color_gradient(low = "steelblue", high = "firebrick") +
  labs(
    title = "KEGG enrichment – Vegetative system",
    x     = "Gene ratio (DEGs / pathway genes)",
    y     = "KEGG pathway",
    size  = "DEG count",
    color = "-log10(adj p-value)"
  ) +
  theme_minimal(base_size = 12)