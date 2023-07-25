
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggtree)
    library(treeio)

})

#setwd("~/Documents/GEOMAR/genotyping/mitotype_pipeline_revision")

grouping <- read_tsv("./labels_all.txt") %>% distinct()

tree <- read.nexus(file="./output/tonsa_mb.nex.con.tre")
x <- as_tibble(tree)

dm <- left_join(x, grouping, by="label")



p <- ggtree(tree) + theme_tree()


pout <- p %<+% grouping + 
    geom_tiplab(aes(color=group), size=0.9) +
    theme(legend.position="right")+ 
    #geom_text( show.legend  = F ) +
    geom_tippoint(aes(color=group), size=0.9) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = c("A_hudsonica" = "#1B9E77",
                               "A_lilljeborgi" = "#D95F02",
                               "F" = "#7570B3",
                               "IV" = "#E7298A",
                               "out_group" = "#666666",
                               "S" = "#66A61E",
                               "SB" = "#E6AB02",
                               "X" = "#A6761D"),
                    na.value = "black")


pout <- p %<+% grouping + 
  geom_tiplab(aes(color=group), size=0.9) +
  theme(legend.position="right")+ 
  #geom_text( show.legend  = F ) +
  geom_tippoint(aes(color=group), size=0.9) +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values = c("A_hudsonica" = "#1B9E77",
                                "A_tonsa" = "#D95F02",
                                "outlier" = "firebrick2",
                                "Collection_1" = "darkorchid2",
                                "Collection_2" = "#66A61E",
                                "Collection_3" = "#E6AB02",
                                "Collection_4" = "#A6761D",
                                "Collection_5" = "dodgerblue2"),
                     na.value = "black")

  



ggsave(pout, file="./output/tree_plot_new.pdf", h=13, w=10)

