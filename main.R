library(Biobase)
library(BiocManager)
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

write_results = FALSE
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
# keep <- rowSums(counts(dds_dson)) > 1
keep <- rowSums(counts(dds_dson)>=5) >= 5
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

# You'll need to have the 'apeglm' package for this to work properly.
res_ae_1 <- generate_deseq_results(dds., name="treatment_AE_1_vs_Baseline", coef=2)
res_ae_4 <- generate_deseq_results(dds., name="treatment_AE_4_vs_Baseline", coef=3)
res_re_1 <- generate_deseq_results(dds., name="treatment_RE_1_vs_Baseline", coef=4)
res_re_4 <- generate_deseq_results(dds., name="treatment_RE_4_vs_Baseline", coef=5)
summary(res_ae_1$apeglm_results, alpha =0.05)
summary(res_ae_4$apeglm_results, alpha =0.05)
summary(res_re_1$apeglm_results, alpha =0.05)
summary(res_re_4$apeglm_results, alpha =0.05)

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
ae.1hr_df <- DE_results(res_ae_1$apeglm_results)
ae.1hr_df_anno <- my_annotations_dickinson(ae.1hr_df)
final.ae_1hr <- left_join(ae.1hr_df, ae.1hr_df_anno,
                          by=c("ensgene"="ensembl_gene_id")) %>%
    dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)

ae.4hr_df <- DE_results(res_ae_4$apeglm_results)
ae.4hr_df_anno <- my_annotations_dickinson(ae.4hr_df)
final.ae_4hr <- left_join(ae.4hr_df, ae.4hr_df_anno,
                          by=c("ensgene"="ensembl_gene_id")) %>%
    dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)

re.1hr_df <- DE_results(res_re_1$apeglm_results)
re.1hr_df_anno <- my_annotations_dickinson(re.1hr_df)
final.re_1hr <- left_join(re.1hr_df, re.1hr_df_anno,
                          by=c("ensgene"="ensembl_gene_id")) %>%
    dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)

re.4hr_df <- DE_results(res_re_4$apeglm_results)
re.4hr_df_anno <- my_annotations_dickinson(re.4hr_df)
final.re_4hr <- left_join(re.4hr_df, re.4hr_df_anno,
                          by=c("ensgene"="ensembl_gene_id")) %>%
    dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)

print_gene_results <- function(){
    cat("---1HR Aerobic Exercise Top DE genes---\n")
    final.ae_1hr %>% head %>% print
    cat("\n---4HR Aerobic Exercise Top DE genes---\n")
    final.ae_4hr %>% head %>% print
    cat("\n---1HR Resistance Exercise Top DE genes---\n")
    final.re_1hr %>% head %>% print
    cat("\n---4HR Resistance Exercise Top DE genes---\n")
    final.re_4hr %>% head %>% print
}
print_gene_results() # from the supplementary

print_lncRNA_tables <- function() {
    cat("---1HR Aerobic Exercise Top DE LncRNAs---\n")
    final.ae_1hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        print
    cat("\n---4HR Aerobic Exercise Top DE LncRNAs---\n")
    final.ae_4hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        print
    cat("\n---1HR Resistance Exercise Top DE LncRNAs---\n")
    final.re_1hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        print
    cat("\n---4HR Resistance Exercise Top DE LncRNAs---\n")
    final.re_4hr %>% filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" ) %>%
        print
}
print_lncRNA_tables() # from the main-text

final.ae_1hr <- final.ae_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)
final.ae_4hr <- final.ae_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)
final.re_1hr <- final.re_1hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)
final.re_4hr <- final.re_4hr %>% dplyr::select(ensgene, external_gene_name, log2FoldChange, lfcSE, pvalue, padj, gene_biotype)

if (!file.exists("proc")) dir.create("proc")
if (write_results == TRUE) {
    write_csv(final.ae_1hr, file = "proc/AE_baseline.vs.1.csv", append = FALSE)
    write_csv(final.ae_4hr, file = "proc/AE_baseline.vs.4.csv", append = FALSE)
    write_csv(final.re_1hr, file = "proc/RE_baseline.vs.1.csv", append = FALSE)
    write_csv(final.re_4hr, file = "proc/RE_baseline.vs.4.csv", append = FALSE)
}

## Extra plots
# Volcano plot

df_ae_1 <- as.data.frame(res_ae_1$normal_results) %>%
    dplyr::mutate(threshold = padj<0.05 & abs(log2FoldChange) >= 0.58) %>%
    subset(!is.na(baseMean)) %>%
    rownames_to_column(var = "ensgene")

df_ae_4 <- as.data.frame(res_ae_4$normal_results) %>%
    dplyr::mutate(threshold = padj<0.05 & abs(log2FoldChange) >= 0.58) %>%
    subset(!is.na(baseMean)) %>%
    rownames_to_column(var = "ensgene")

df_re_1 <- as.data.frame(res_re_1$normal_results) %>%
    dplyr::mutate(threshold = padj<0.05 & abs(log2FoldChange) >= 0.58) %>%
    subset(!is.na(baseMean)) %>%
    rownames_to_column(var = "ensgene")

df_re_4 <- as.data.frame(res_re_4$normal_results) %>%
    dplyr::mutate(threshold = padj<0.05 & abs(log2FoldChange) >= 0.58) %>%
    subset(!is.na(baseMean)) %>%
    rownames_to_column(var = "ensgene")

anno_DF<- function (df) {
    anno <- my_annotations_dickinson(tdf)
    anno_gene <- dplyr::distinct(anno)
    df$anno <- anno[match(df$ensgene, anno_gene$ensembl_gene_id),]$external_gene_name
    df$biotype <- anno[match(df$ensgene, anno_gene$ensembl_gene_id),]$gene_biotype
    df$annoplot <- ifelse(df$threshold == TRUE & (df$biotype == "lincRNA" | df$biotype == "antisense" | df$biotype == "sense_intronic"),
                          df$anno,
                           "")
    df
}

df_ae_1 <- anno_DF(df_ae_1)
df_ae_4 <- anno_DF(df_ae_4)
df_re_1 <- anno_DF(df_re_1)
df_re_4 <- anno_DF(df_re_4)


#anno[tdf$ensgene %in% anno_gene$ensembl_gene_id,]$external_gene_name
# annoplot <- ifelse(tdf$threshold == FALSE, "", tdf$anno)

p1 <- ggplot(df_ae_1,aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color=threshold)) +
    geom_text_repel(aes(label = annoplot), max.overlaps = 100,
                    point.size=NA,
                    segment.colour = NA) +
    ggtitle("1hr") +
    xlim(c(-10, 10)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))

p2 <- ggplot(df_ae_4,aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color=threshold)) +
    geom_text_repel(aes(label = annoplot), max.overlaps = 100,
                    point.size=NA,
                    segment.colour = NA) +
    ggtitle("4hr") +
    xlim(c(-10, 10)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))

(p1 | p2) + plot_annotation(
    title = 'Long Non-Coding RNAs after Aerobic Exercise',
)

p3 <- ggplot(df_re_1,aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color=threshold)) +
    geom_text_repel(aes(label = annoplot), max.overlaps = 100,
                    point.size=NA,
                    segment.colour = NA) +
    ggtitle("1hr") +
    xlim(c(-10, 10)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))

p4 <- ggplot(df_re_4,aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color=threshold)) +
    geom_text_repel(aes(label = annoplot), max.overlaps = 100,
                    point.size=NA,
                    segment.colour = NA) +
    ggtitle("4hr") +
    xlim(c(-10, 10)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))

(p3 | p4) + plot_annotation(
    title = 'Long Non-Coding RNAs after Resistance Exercise',
)

if (!file.exists("images")) dir.create(images)

# Heatmap
