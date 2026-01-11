install.packages(c("GenomicRanges", "data.table"))

library(GenomicRanges)
library(data.table)

shared_deg_bed <- fread(
  "Results/Study_Comparison/Cross_Stress_Comparison/shared_DEGs_B73v5.bed",
  col.names = c("chr", "start", "end", "gene_id")
)

#Sanity Check 
nrow(shared_deg_bed)
head(shared_deg_bed)

background_bed <- fread(
  "B73v5_genes.bed",
  col.names = c("chr", "start", "end", "gene_id")
)

#Loading GWAS Data 

yield_gwas <- fread(
  "GWAS/NAM_GWAS_yield_50kb.bed",
  col.names = c("chr", "start", "end", "trait")
)

flowering_gwas <- fread(
  "GWAS/NAM_GWAS_flowering_50kb.bed",
  col.names = c("chr", "start", "end", "trait")
)

#Convert everything to GenomicRanges
gr_shared <- GRanges(
  seqnames = shared_deg_bed$chr,
  ranges = IRanges(shared_deg_bed$start + 1, shared_deg_bed$end),
  gene_id = shared_deg_bed$gene_id
)

gr_background <- GRanges(
  seqnames = background_bed$chr,
  ranges = IRanges(background_bed$start + 1, background_bed$end),
  gene_id = background_bed$gene_id
)

gr_yield <- GRanges(
  seqnames = yield_gwas$chr,
  ranges = IRanges(yield_gwas$start + 1, yield_gwas$end)
)

gr_flowering <- GRanges(
  seqnames = flowering_gwas$chr,
  ranges = IRanges(flowering_gwas$start + 1, flowering_gwas$end)
)

#Shared DEGs overlapping GWAS loci
shared_yield_hits <- unique(
  mcols(gr_shared)$gene_id[
    queryHits(findOverlaps(gr_shared, gr_yield))
  ]
)

shared_flowering_hits <- unique(
  mcols(gr_shared)$gene_id[
    queryHits(findOverlaps(gr_shared, gr_flowering))
  ]
)

length(shared_yield_hits)
length(shared_flowering_hits)


background_yield_hits <- unique(
  mcols(gr_background)$gene_id[
    queryHits(findOverlaps(gr_background, gr_yield))
  ]
)

background_flowering_hits <- unique(
  mcols(gr_background)$gene_id[
    queryHits(findOverlaps(gr_background, gr_flowering))
  ]
)

#Enrichment tests 
yield_table <- matrix(
  c(
    length(shared_yield_hits),
    length(gr_shared) - length(shared_yield_hits),
    length(background_yield_hits),
    length(gr_background) - length(background_yield_hits)
  ),
  nrow = 2,
  byrow = TRUE
)

colnames(yield_table) <- c("GWAS_overlap", "No_overlap")
rownames(yield_table) <- c("Shared_DEGs", "Background")

yield_fisher <- fisher.test(yield_table)

yield_table
yield_fisher

flowering_table <- matrix(
  c(
    length(shared_flowering_hits),
    length(gr_shared) - length(shared_flowering_hits),
    length(background_flowering_hits),
    length(gr_background) - length(background_flowering_hits)
  ),
  nrow = 2,
  byrow = TRUE
)

colnames(flowering_table) <- c("GWAS_overlap", "No_overlap")
rownames(flowering_table) <- c("Shared_DEGs", "Background")

flowering_fisher <- fisher.test(flowering_table)

flowering_table
flowering_fisher

saveRDS(shared_yield_hits, "Results/shared_DEGs_yieldGWAS_hits.rds")
saveRDS(shared_flowering_hits, "Results/shared_DEGs_floweringGWAS_hits.rds")

#Visualize using BArplot 
library(ggplot2)

enrichment_df <- data.frame(
  Category = c("Yield GWAS", "Yield GWAS", "Flowering GWAS", "Flowering GWAS"),
  Group = c("Shared DEGs", "Background", "Shared DEGs", "Background"),
  Proportion = c(
    65/165,
    12696/40140,
    28/165,
    7645/40140
  )
)

ggplot(enrichment_df, aes(x = Category, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    y = "Proportion of genes overlapping GWAS loci",
    x = "",
    title = "Enrichment of shared drought–cold DEGs near GWAS loci"
  ) +
  theme_minimal(base_size = 13)


library(GenomicRanges)

# Label which shared DEGs overlap yield GWAS
shared_deg_bed$Yield_GWAS <- shared_deg_bed$gene_id %in% shared_yield_hits

ggplot(shared_deg_bed, aes(
  x = start / 1e6,
  y = chr,
  color = Yield_GWAS
)) +
  geom_point(alpha = 0.7) +
  labs(
    x = "Genomic position (Mb)",
    y = "Chromosome",
    title = "Genomic distribution of shared drought–cold DEGs"
  ) +
  scale_color_manual(
    values = c("grey70", "red"),
    labels = c("No yield GWAS overlap", "Yield GWAS overlap")
  ) +
  theme_minimal(base_size = 12)
