library(ggvenn)
library(RColorBrewer)

write_image = FALSE

get_lnc_names <- function(df) {
    df %>%
        filter(gene_biotype == "lincRNA" | gene_biotype == "antisense" | gene_biotype == "sense_intronic") %>%
        dplyr::select(external_gene_name) %>%
        unlist
}

ae_1_4 <- list(`1 Hour`=ae_1_lnc,
               `4 Hours`=ae_4_lnc)

re_1_4 <- list(`1 Hour`=re_1_lnc,
                `4 Hours`=re_4_lnc)

all_ex <- list(`1 Hour RE`=re_1_lnc,
               `4 Hours RE`=re_4_lnc,
               `1 Hour AE`=ae_1_lnc,
               `4 Hours AE`=ae_4_lnc)

ae_re_long <- list(`Aerobic`=ae_4_lnc,
              `Resistance`=re_4_lnc)

ae_re_short <- list(`Aerobic`=ae_1_lnc,
                   `Resistance`=re_1_lnc)


make_venn <- function(listdata, title, fillname, saveImage, showElements=TRUE, nCol=3) {
    ggimage <- listdata %>%
        ggvenn(show_elements = showElements,
               label_sep = "\n",
               fill_color = brewer.pal(name=fillname,n=nCol)) +
        ggtitle(title) +
        theme(plot.title = element_text(color = "#0099f8", size = 18, face = "bold"))
    print(ggimage)
    if (saveImage == TRUE) ggsave(paste0("images/",gsub(" ", "_", title), ".png"))
}
make_venn(all_ex, title="Exercise Modality", fillname="Set3", saveImage=write_image, showElements=FALSE, nCol=4)
make_venn(ae_1_4, title="Aerobic Exercise", fillname="Set2", saveImage=write_image, showElements=TRUE, nCol=3)
make_venn(re_1_4, title="Resistance Exercise", fillname="Set2", saveImage=write_image, showElements=TRUE, nCol=3)
make_venn(ae_re_short, title="1 Hour Exercise", fillname="Set2", saveImage=write_image, showElements=TRUE, nCol=3)
make_venn(ae_re_long, title="4 Hour Exercise", fillname="Set2", saveImage=write_image, showElements=TRUE, nCol=3)
