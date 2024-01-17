library(Biobase)
library(BiocManager)
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)

pdata_dson <- getGEO("GSE107934")[[1]]
citation("GEOquery")

pheno_dson <- pData(pdata_dson) %>%
    dplyr::select(geo_accession,
                  `exercise:ch1`,
                  `time:ch1`)
pheno_dson_clean <- dplyr::rename(pheno_dson,
                                  all_of(
                                        c(exercise="exercise:ch1", time="time:ch1")
                                      )
                                  )
## Annotating conditions and treatments
treatment_dson <- c("Baseline", "AE_1", "AE_4", "RE_1","RE_4")
rep_treat_dson <- rep(treatment_dson,6)[-18] # removes the 18th entry ("AE_4") as it was not in the final analysis
factor_treat <- relevel(factor(rep_treat_dson), ref = "Baseline")
pheno_dson_clean$treatment <- factor_treat
pheno_dson_clean

## Getting count data
if (!file.exists("raw/GSE107934_RAW.tar")) {
    dir.create("raw")
    getGEOSuppFiles("GSE107934", baseDir = "raw", makeDirectory = FALSE)
}
if (!file.exists("raw/dson")) untar("raw/GSE107934_RAW.tar", exdir = "raw/dson")

my_files <- list.files("raw/dson", full.names=TRUE)
count_read <- lapply(my_files, read.table)
count_df <- as.data.frame(count_read)

# Formatting dataframe
rownames(count_df) <- count_df$V1
new_count_df <- count_df %>%
    dplyr::select(everything()[c(FALSE,TRUE)])

cts_dson <- as.matrix(new_count_df)
cts_dson
colnames(cts_dson) = rownames(pheno_dson_clean)
colnames(cts_dson)

## Creating DESEQDataSet
dds_dson <- DESeqDataSetFromMatrix(countData = cts_dson,
                       colData = pheno_dson_clean,
                       design = ~treatment)

# Filtering all genes that have less than 5 reads across samples
keep <- rowSums(counts(dds_dson)) > 1
# keep <- rowSums(counts(dds_dson) >= 5)
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

generate_deseq_results <- function(dds., name=name, coef=coef) {
    res_norm <- results(dds., name=name)
    shrunk_res <- lfcShrink(dds., coef=coef, res=res_norm)
    return(list(normal_results=res_norm, apeglm_results=shrunk_res, exp_name=name))
}

generate_MA_plots <- function(res, xlim=c(1,1e5), ylim=c(-11,11)) {
    par(mfrow=c(1,2))
    plotMA(res$normal_results, ylim=ylim, main="Normal")
    plotMA(res$apeglm_results,ylim=ylim, main="LFC apeglm")
    mtext(res$exp_name, side = 3, line = -1, outer = TRUE, cex=1,font=2)
}

res_ae_1 <- generate_deseq_results(dds., name="treatment_AE_1_vs_Baseline", coef=2)
res_ae_4 <- generate_deseq_results(dds., name="treatment_AE_4_vs_Baseline", coef=3)
res_re_1 <- generate_deseq_results(dds., name="treatment_RE_1_vs_Baseline", coef=4)
res_re_4 <- generate_deseq_results(dds., name="treatment_RE_4_vs_Baseline", coef=5)
summary(res_ae_1$apeglm_results, alpha =0.05)

# AE 1 hour
generate_MA_plots(res=res_ae_1)
generate_MA_plots(res=res_ae_4)
generate_MA_plots(res=res_re_1)
generate_MA_plots(res=res_re_4)

DE_results = function(DESeqResults_class) {
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
my_ae.1hr_df <- DE_results(res_ae_1$apeglm_results)
my_ae.1hr_df_anno <- my_annotations_dickinson(my_ae.1hr_df)
final.ae_1hr <- left_join(my_ae.1hr_df, my_ae.1hr_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

ae_1.4_my.r_df <- DE_results(res_ae_4$apeglm_results)
ae_1.4_my.r_df_anno <- my_annotations_dickinson(ae_1.4_my.r_df)
final.ae_4hr <- left_join(ae_1.4_my.r_df, ae_1.4_my.r_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

re_b.1.r_df <- DE_results(res_re_1$apeglm_results)
re_b.1.r_df_anno <- my_annotations_dickinson(re_b.1.r_df)
final.re_1hr <- left_join(re_b.1.r_df, re_b.1.r_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

re_1.4.r_df <- DE_results(res_re_4$apeglm_results)
re_1.4.r_df_anno <- my_annotations_dickinson(re_1.4.r_df)
final.re_4hr <- left_join(re_1.4.r_df, re_1.4.r_df_anno,
                          by=c("ensgene"="ensembl_gene_id"))

print_gene_results <- function(){
    cat("---1HR Aerobic Exercise Top DE genes---\n")
    final.ae_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
    cat("\n---4HR Aerobic Exercise Top DE genes---\n")
    final.ae_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
    cat("\n---1HR Resistance Exercise Top DE genes---\n")
    final.re_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
    cat("\n---4HR Resistance Exercise Top DE genes---\n")
    final.re_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head %>% print
}
print_gene_results() # from the supplementary

print_lncRNA_tables <- function() {
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
print_lncRNA_tables() # from the main-text

final.ae_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head
final.ae_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head
final.re_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head
final.re_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype) %>% head

