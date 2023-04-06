library(reshape2)
library(vegan)
library(ecodist)
library(ape)
library(plyr)
library(dplyr)
library(colortools)
library(data.table)
library(factoextra)
library(ggplot2)
library(corrplot)
library(vegan)

setwd("/path/to/github/repo/NICUBSI/tables")

sample_df = read.delim("shortbred_mg.tsv")
metadata = read.csv("MasterMetadata_FINAL_GitHub.csv")

sample_df = rbind(sample_cat, sample_not_cat)
sample_df$exp_type = metadata$exp_type[match(sample_df$sample, metadata$Sequence.name)]
sample_df = sample_df[sample_df$exp_type %in% c("control", "experimental"),]
sample_df$dpi = metadata$DPI1[match(sample_df$sample, metadata$Sequence.name)]
sample_df$subject = metadata$Subject[match(sample_df$sample, metadata$Sequence.name)]
sample_df$before_infx = sample_df$dpi < 0
sample_df$before_infx[sample_df$exp_type == "control"] = T
sample_df$group = paste(sample_df$subject, sample_df$before_infx)

cast_mg_abx = dcast(sample_df, sample ~ gene, value.var="relative_abundance", fun.aggregate=length)
cast_mg_abx[is.na(cast_mg_abx)] = 0

rownames(cast_mg_abx) = cast_mg_abx$sample

cast_mg_abx = cast_mg_abx[, 2:ncol(cast_mg_abx)]
cast_mg_abx$exp_type = metadata$exp_type[match(rownames(cast_mg_abx), metadata$Sequence.name)]
cast_mg_abx$dpi = metadata$DPI1[match(rownames(cast_mg_abx), metadata$Sequence.name)]
cast_mg_abx$subject = metadata$Subject[match(rownames(cast_mg_abx), metadata$Sequence.name)]
cast_mg_abx$before_infx = cast_mg_abx$dpi < 0
cast_mg_abx$group = paste(cast_mg_abx$subject, cast_mg_abx$exp_type)
cast_mg_abx = cast_mg_abx[!is.na(cast_mg_abx$exp_type),]

df_list = lapply(unique(as.character(cast_mg_abx$group)), function(g) {
  df = cast_mg_abx[cast_mg_abx$group == g, 1:(ncol(cast_mg_abx)-4)]
  mean_df = colSums(df) / dim(df)[1]
  as.list(mean_df)
})

df_aggregated = as.data.frame(rbindlist(df_list))

d = df_aggregated[!is.na(cast_mg_abx$exp_type[!duplicated(cast_mg_abx$group)]), as.numeric(sapply(df_aggregated, var)) > 0]

pca.abx.genes = prcomp(d, scale = TRUE)

fviz_pca_ind(
    pca.abx.genes,
    col.ind = experimental_groups, # color by groups
    geom = "point",
    title = "Subject avg shortbred ARGs",
    pointshape = "circle",
    palette = c("grey50", "red", "blue"),
    addEllipses = TRUE, # Concentration ellipses
    ellipse.type = "confidence",
    legend.title = "Groups",
    mean.point = TRUE,
    repel = TRUE
)
