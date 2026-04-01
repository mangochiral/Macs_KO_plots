library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(ggrepel)

res_A20_vs_AAVS1

# Create dataframe with gene symbols
# Create dataframe - no need for rownames_to_column or left_join
de_results <- res_A20_vs_AAVS1 %>%
  as.data.frame() %>%
  filter(!is.na(padj))

# Get significant genes (for enrichment)
sig_genes <- de_results %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(gene_symbol)

cat("Number of significant genes:", length(sig_genes), "\n")


# 2) Functional Enrichment Analysis - GO
# ==============================================================================

# GO Biological Process enrichment
ego_BP <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# View top enriched pathways
head(ego_BP, 20)

# Save GO results
ego_BP %>%
  as.data.frame() %>%
  write_csv("output/GO_BP_enrichment_results.csv")

# Plot GO enrichment
dotplot(ego_BP, showCategory = 15, title = "GO Biological Process Enrichment") +
  theme(axis.text.y = element_text(size = 10))
ggsave("plots/02.01_GO_BP_dotplot.png", width = 10, height = 8)

# ==============================================================================
# 3) KEGG Pathway Enrichment
# ==============================================================================

# Convert gene symbols to Entrez IDs for KEGG
gene_entrez <- bitr(sig_genes, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# KEGG enrichment
ekegg <- enrichKEGG(
  gene         = gene_entrez$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

head(ekegg, 20)

# Save KEGG results
ekegg %>%
  as.data.frame() %>%
  write_csv("output/KEGG_enrichment_results.csv")

# Plot KEGG enrichment
dotplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment") +
  theme(axis.text.y = element_text(size = 10))
ggsave("plots/02.02_KEGG_dotplot.png", width = 10, height = 8)

# ==============================================================================
# 4) MSigDB Hallmark Gene Sets
# ==============================================================================
# Get Hallmark gene sets - use dplyr::select explicitly
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# Enrichment with Hallmark gene sets
ehallmark <- enricher(
  gene = sig_genes,
  TERM2GENE = hallmark_sets,
  pvalueCutoff = 0.05
)

head(ehallmark, 20)

# Save Hallmark results
ehallmark %>%
  as.data.frame() %>%
  write_csv("output/Hallmark_enrichment_results.csv")

# Plot Hallmark enrichment
dotplot(ehallmark, showCategory = 15, title = "Hallmark Gene Sets Enrichment") +
  theme(axis.text.y = element_text(size = 10))
ggsave("plots/02.03_Hallmark_dotplot.png", width = 10, height = 8)
