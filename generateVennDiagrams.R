library(VennDiagram)
get_lnc_names <- function(df) {
    df %>%
        filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        dplyr::select(external_gene_name) %>%
        unlist
}

ae_1_lnc <- get_lnc_names(final.ae_1hr)
ae_4_lnc <- get_lnc_names(final.ae_4hr)
re_1_lnc <- get_lnc_names(final.re_1hr)
re_4_lnc <- get_lnc_names(final.re_4hr)

ae_lnc <- union(ae_1_lnc, ae_4_lnc)
re_lnc <- union(re_1_lnc, re_4_lnc)

get_sets <- function(set1, set2) {
    total <- union(set1, set2)
    cross <- length(intersect(set1, set2))
    right <- length(setdiff(total, set1))
    left <- length(total)-right
    return(list(left=left, right=right, cross=cross))
}

ae_1v4 <- get_sets(ae_1_lnc, ae_4_lnc)
re_1v4 <- get_sets(re_1_lnc, re_4_lnc)
ae_v_re <- get_sets(ae_lnc, re_lnc)

grid.newpage()
venn.plot <- draw.pairwise.venn(ae_1v4$left, ae_1v4$right, ae_1v4$cross, category=c("1hr Aerobic", "4hr Aerobic"))
grid.draw(venn.plot)

grid.newpage()
venn.plot <- draw.pairwise.venn(re_1v4$left, re_1v4$right, re_1v4$cross, category=c("1hr Resistance", "4hr Resistance"))
grid.draw(venn.plot)

grid.newpage()
venn.plot <- draw.pairwise.venn(ae_v_re$left, ae_v_re$right, ae_v_re$cross, category=c("Aerobic", "Resistance"))
grid.draw(venn.plot)
