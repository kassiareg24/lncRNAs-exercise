make_ego <- function (res, tx2gene) {
  genes <- unique(tx2gene[tx2gene$gene_name %in% res$gene_name,]$gene_id) # gets unique genes from symbols
  ego <- enrichGO(gene = genes,
                  universe = tx2gene$gene_id,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'ENSEMBL',
                  readable = T,
                  ont = "BP",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.10)
}

library(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
unique_RNA = unique(c(final.ae_1hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>% pull(ensgene), 
  final.ae_4hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>% pull(ensgene),
  final.re_1hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>% pull(ensgene), 
  final.re_4hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>% pull(ensgene)))

mappable <- unique_RNA %in% keys(org.Hs.eg.db, keytype = 'ENSEMBL')
mappedLNC <- unique_RNA[mappable]
DE_AE <- c(final.ae_1hr$ensgene, final.ae_4hr$ensgene)
DE_RE <- c(final.re_1hr$ensgene, final.re_4hr$ensgene)

GO_AE <- enrichGO(gene = DE_AE,
         universe = rownames(dds),
         OrgDb = "org.Hs.eg.db",
         keyType = 'ENSEMBL',
         readable = T,
         ont = "BP",
         pAdjustMethod="BH")

GO_RE <- enrichGO(gene = DE_RE,
         universe = rownames(dds),
         OrgDb = "org.Hs.eg.db",
         keyType = 'ENSEMBL',
         readable = T,
         ont = "BP",
         pAdjustMethod="BH")

GO_AE_simple <- simplify(GO_AE, cutoff=0.7, by="p.adjust")
GO_RE_simple <- simplify(GO_RE, cutoff=0.7, by="p.adjust")

debug(dotplot)
dotplot(GO_AE_simple, showCategory=15, title="Biological Processes for Aerobic Exercise", color='p.adjust')
undebug(dotplot)
dotplot(GO_RE_simple, showCategory=15, title="Biological Processes for Resistance Exercise", color='p.adjust')
## ----

ae_df <- GO_AE_simple@result[, c("Description", "GeneRatio", "geneID")]
re_df <- GO_RE_simple@result[, c("Description", "GeneRatio", "geneID")]

specific_term_ae <- ae_df[ae_df$Description == "positive regulation of cytokine production", ]
specific_term_re <- re_df[re_df$Description == "positive regulation of cytokine production", ]

specific_term_re
specific_term_ae

specific_term_re <- re_df[re_df$Description == "myeloid leukocyte differentiation", ]
specific_term_ae <- ae_df[ae_df$Description == "myeloid leukocyte differentiation", ]

specific_term_re
specific_term_ae

specific_term_re = re_df[re_df$Description == "skeletal muscle cell differentiation", ]
specific_term_ae = ae_df[ae_df$Description == "skeletal muscle cell differentiation", ]

specific_term_re
specific_term_ae
