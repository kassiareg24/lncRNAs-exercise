library(Biobase)
library(BiocManager)
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)

pdata_dson <- getGEO("GSE107934")
pdata_dson <- pdata_dson[[1]]
citation("GEOquery")

pheno_dson <- pData(pdata_dson) %>%
    dplyr::select(geo_accession,
                  `exercise:ch1`,
                  `time:ch1`)

pheno_dson$exercise <- pheno_dson$`exercise:ch1`
pheno_dson$time <- pheno_dson$`time:ch1`

pheno_dson_clean = subset(pheno_dson, select= -c(`exercise:ch1`, `time:ch1`))

treatment_dson <- c("Baseline", "AE_1", "AE_4", "RE_1","RE_4")
rep_treat_dson <- c(rep("Baseline",1), rep("AE_1",1), rep("AE_4",1), rep("RE_1",1),rep("RE_4",1),
                    rep("Baseline",1), rep("AE_1",1), rep("AE_4",1), rep("RE_1",1),rep("RE_4",1),
                    rep("Baseline",1), rep("AE_1",1), rep("AE_4",1), rep("RE_1",1),rep("RE_4",1),
                    rep("Baseline",1), rep("AE_1",1), rep("RE_1",1), rep("RE_4",1),
                    rep("Baseline",1), rep("AE_1",1), rep("AE_4",1), rep("RE_1",1),rep("RE_4",1),
                    rep("Baseline",1), rep("AE_1",1), rep("AE_4",1), rep("RE_1",1),rep("RE_4",1))

factor_treat <- factor(rep_treat_dson)
pheno_dson_clean = subset(pheno_dson, select = -c(`exercise:ch1`,`time:ch1`))

pheno_dson_clean$treatment <- factor_treat
releveld_treat <- relevel(pheno_dson_clean$treatment, ref = "Baseline")
pheno_dson_clean$treatment = releveld_treat
getGEOSuppFiles("GSE107934")
untar("GSE107934/GSE107934_RAW.tar", exdir = "dson")

my_files <- list.files("dson", full.names=TRUE)
count_read <- lapply(my_files, read.table)
count_df <- as.data.frame(count_read)
rownames(count_df) <- count_df$V1

new_count_df <- count_df %>%
    dplyr::select(everything()[c(FALSE,TRUE)])

# new_count_df <- new_count_df[-seq(2,58,by=2)] # is this line necessary?

length(seq(2,58,by=2))
cts_dson <- as.matrix(new_count_df)
cts_dson
colnames(cts_dson) = rownames(pheno_dson_clean)
colnames(cts_dson)

## Creating DESEQDataSet
dds_dson <- DESeqDataSetFromMatrix(countData = cts_dson,
                       colData = pheno_dson_clean,
                       design = ~treatment)

keep <- rowSums(counts(dds_dson)) > 1
dds <-dds_dson[keep,]

dds. <- DESeq(dds)

resultsNames(dds.)
boxplot(counts(dds., normalized = TRUE))

vst <- varianceStabilizingTransformation(dds.)
boxplot(assay(vst))

plotPCA(vst, intgroup = c("time", "exercise"))

pcaData <- plotPCA(vst, intgroup = c("time", "exercise"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x=PC1, y=PC2, color=time, shape=exercise)) +
    geom_point(size=3) +
    xlab(paste0("PC1:   ", percentVar[1], "% variance")) +
    ylab(paste0("PC2:   ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle("PCA with VST data")

res.ae.b.1 <- results(dds., name="treatment_AE_1_vs_Baseline")
summary(res.ae.b.1)

res.ae.b.4 <- results(dds., name = "treatment_AE_4_vs_Baseline")
res.re.b.1 <- results(dds., name = "treatment_RE_1_vs_Baseline")
res.re.b.4 <- results(dds., name = "treatment_RE_4_vs_Baseline")

resLFC_ae.1 <- lfcShrink(dds.,
                           coef=2,
                           res=res.ae.b.1)

resLFC_ae.b.4 <- lfcShrink(dds.,
                           coef=3,
                           res=res.ae.b.4)

resLFC_re.b.1 <- lfcShrink(dds.,
                           coef=4,
                           res=res.re.b.1)

resLFC_re.1.4 <- lfcShrink(dds.,
                           coef=5,
                           res=res.re.b.4)

xlim <- c(1, 1e5)
ylim <- c(-11,11)

# AE 1 hour
plotMA(res.ae.b.1, ylim=ylim, main="Normal")
plotMA(resLFC_ae.1,ylim=ylim, main="apeglm")

# AE 4 hours
plotMA(res.ae.b.4,ylim=ylim,main="Normal")
plotMA(resLFC_ae.b.4,ylim=ylim,main="apeglm")

# RE 1 hour
plotMA(res.re.b.1,ylim=ylim,main="Normal")
plotMA(resLFC_re.b.1,xlim=xlim,ylim=ylim,main="apeglm")

# RE 4 hours
plotMA(res.re.b.4,ylim=ylim,main="Normal")
plotMA(resLFC_re.1.4,ylim=ylim,main="apeglm")

my_results = function(DESeqResults_class) {
    df = as.data.frame(DESeqResults_class)
    filter_padj = df[df$padj < 0.05,]
    filter_logFC = filter_padj[abs(filter_padj$log2FoldChange) > 0.58,]
    clean_df = subset(filter_logFC, !is.na(baseMean))
    r2c=rownames_to_column(clean_df, var = "ensgene")
    return(r2c)
}
my_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   host="https://grch37.ensembl.org",
                   path="/biomart/martservice",
                   dataset="hsapiens_gene_ensembl")

my_annotations_dickinson = function(x3) {
    y= getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name",
                            "gene_biotype",
                            "chromosome_name",
                            "start_position",
                            "end_position",
                            "strand",
                            "description"),
             filters= c("ensembl_gene_id"), values=x3$ensgene,
             mart = my_mart)
}

# Note, if this fails, use the following version of dbplyr:
# devtools::install_version("dbplyr", version = "2.3.4")

resLFC_ae.1hr_df <- as.data.frame(resLFC_ae.1)
my_ae.1hr_df <- my_results(resLFC_ae.1hr_df)
my_ae.1hr_df_anno <- my_annotations_dickinson(my_ae.1hr_df)
final.ae_1hr <- left_join(my_ae.1hr_df, my_ae.1hr_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

ae_1.4_my.r_df <- my_results(resLFC_ae.b.4)
ae_1.4_my.r_df_anno <- my_annotations_dickinson(ae_1.4_my.r_df)
final.ae_4hr <- left_join(ae_1.4_my.r_df, ae_1.4_my.r_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

re_b.1.r_df <- my_results(resLFC_re.b.1)
re_b.1.r_df_anno <- my_annotations_dickinson(re_b.1.r_df)
final.re_1hr <- left_join(re_b.1.r_df, re_b.1.r_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

re_1.4.r_df <- my_results(resLFC_re.1.4)
re_1.4.r_df_anno <- my_annotations_dickinson(re_1.4.r_df)
final.re_4hr <- left_join(re_1.4.r_df, re_1.4.r_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

print_results <- function(){
    cat("---1HR Aerobic Exercise Top DE genes---\n")
    final.ae_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
    cat("\n---4HR Aerobic Exercise Top DE genes---\n")
    final.ae_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
    cat("\n---1HR Resistance Exercise Top DE genes---\n")
    final.re_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
    cat("\n---4HR Resistance Exercise Top DE genes---\n")
    final.re_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
}
print_results()

print_tables <- function() {
    cat("---1HR Aerobic Exercise Top DE LncRNAs---\n")
    final.ae_1hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% print
    cat("\n---4HR Aerobic Exercise Top DE LncRNAs---\n")
    final.ae_4hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% print
    cat("\n---1HR Resistance Exercise Top DE LncRNAs---\n")
    final.re_1hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% print
    cat("\n---4HR Resistance Exercise Top DE LncRNAs---\n")
    final.re_4hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" ) %>%
        dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% print
}
print_tables()

final.ae_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head
final.ae_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head
final.re_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head
final.re_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head
