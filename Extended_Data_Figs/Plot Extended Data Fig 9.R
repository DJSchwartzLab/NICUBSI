library(reshape2)
library(vegan)
library(ecodist)
library(ape)
library(plyr)
library(dplyr)
library(ggplot2)
library(colortools)

# load metadata

setwd("/Users/nitan/Documents/Dantas/NeonateBSI/master\ metadata")
metadata = read.csv("MasterMetadata_FINAL.csv")

# load the species relative abundance

# sample specific factors: Abx_days, ABx_score_sample, DOL, NPO_days, HM_days, Form_days,	DM_days, DPI1
# subject specific factors: exp_type, gest_age, birthweight, sex, bug1_Family

species_df = read.delim("/Users/nitan/Documents/Dantas/NeonateBSI/metaphlan/211210_final/211210_wide_metaphlanStandard_species_100k.txt", sep=" ")

# add all the following metadata attributes to the species table
species_df$sample_id = sapply(strsplit(as.character(species_df$Sample), '_'), function(s) s[[1]][1])
species_df$subject = sapply(strsplit(as.character(species_df$Sample), '_'), function(s) s[[1]][1])

species_df$dol = metadata$DOL[match(species_df$Sample, metadata$Sequence.name)]
original_sample_vector = species_df$Sample
species_df = species_df[species_df$dol <= 90,]
removed_samples = original_sample_vector[!(original_sample_vector %in% species_df$Sample)] # samples with DOL > 90

species_df$exp_type = metadata$exp_type[match(species_df$Sample, metadata$Sequence.name)]
species_df$group = metadata$exp_id[match(species_df$Sample, metadata$Sequence.name)]
species_df$ga = metadata$gest_age[match(species_df$Sample, metadata$Sequence.name)]
species_df$bw = metadata$birthweight[match(species_df$Sample, metadata$Sequence.name)]
species_df$dpi = metadata$DPI1[match(species_df$Sample, metadata$Sequence.name)]
species_df$abx_score_complete = metadata$ABx_score_sample[match(species_df$Sample, metadata$Sequence.name)]
species_df$Bug1_Family = metadata$Bug1_Family[match(species_df$Sample, metadata$Sequence.name)]

species_df$Abx_days = metadata$Abx_days[match(species_df$Sample, metadata$Sequence.name)]
species_df$Abx_score_sample = metadata$ABx_score_sample[match(species_df$Sample, metadata$Sequence.name)]
species_df$NPO_days = metadata$NPO_days[match(species_df$Sample, metadata$Sequence.name)]
species_df$Form_days = metadata$Form_days[match(species_df$Sample, metadata$Sequence.name)]
species_df$sex = metadata$sex[match(species_df$Sample, metadata$Sequence.name)]
species_df$DM_days = metadata$DM_days[match(species_df$Sample, metadata$Sequence.name)]
species_df$DPI1 = metadata$DPI1[match(species_df$Sample, metadata$Sequence.name)]
species_df$HM_days = metadata$HM_days[match(species_df$Sample, metadata$Sequence.name)]

# turn rownames to the sample name
rownames(species_df) = species_df$Sample

# load species table again for calculating a distance matrix
species_mat = read.delim("/Users/nitan/Documents/Dantas/NeonateBSI/metaphlan/211210_final/211210_wide_metaphlanStandard_species_100k.txt", sep=" ")
species_mat = species_mat[!(species_mat$Sample %in% removed_samples),]
rownames(species_mat) = species_mat$Sample
# drop samples and species count columns
species = species_mat[, 2:(ncol(species_mat)-1),]

# keep sample names for which there are species counts > 0; i.e. non-empty samples
nonempty_samples = rownames(species)[colSums(t(species)) > 0]

# calculate distance matrix
dis = vegdist(sapply(species, as.numeric),  method = "bray", diagonal=T, na.rm=T)
# calculate PCOA
pcoa_res = pcoa(D=dis)
# extract PCOA dim reduction vectors
df = as.data.frame(pcoa_res$vectors)
pc1 = paste0(as.character(100 * round(pcoa_res$values$Relative_eig[1], 3)), "%")
pc2 = paste0(as.character(100 * round(pcoa_res$values$Relative_eig[2], 3)), "%")

# add the various metadata 
df$exp_type = species_df$exp_type[rownames(species_df) %in% nonempty_samples]
df$sample_name = species_df$Sample[rownames(species_df) %in% nonempty_samples]
df$exp_type = species_df$exp_type[rownames(species_df) %in% nonempty_samples]
df$group = species_df$group[rownames(species_df) %in% nonempty_samples]
df$sample_id = species_df$sample_id[rownames(species_df) %in% nonempty_samples]
df$dol = species_df$dol[rownames(species_df) %in% nonempty_samples]
df$ga = species_df$ga[rownames(species_df) %in% nonempty_samples]
df$bw = species_df$bw[rownames(species_df) %in% nonempty_samples]
df$dpi = species_df$dpi[rownames(species_df) %in% nonempty_samples]
df$dpi = species_df$dpi[rownames(species_df) %in% nonempty_samples]
df$abx_days = species_df$Abx_days[rownames(species_df) %in% nonempty_samples]
df$NPO_days = species_df$NPO_days[rownames(species_df) %in% nonempty_samples]
df$HM_days = species_df$HM_days[rownames(species_df) %in% nonempty_samples]
df$Form_days = species_df$Form_days[rownames(species_df) %in% nonempty_samples]
df$DM_days = species_df$DM_days[rownames(species_df) %in% nonempty_samples]
df$sex = species_df$sex[rownames(species_df) %in% nonempty_samples]
df$sample = rownames(species_df)[rownames(species_df) %in% nonempty_samples]
df$abx_score = species_df$abx_score_complete[rownames(species_df) %in% nonempty_samples]
df$Bug1_Family = species_df$Bug1_Family[rownames(species_df) %in% nonempty_samples]
df$subject = species_df$subject[rownames(species_df) %in% nonempty_samples]
df$class = as.character(df$exp_type)
df$family = df$Bug1_Family

df$class[df$exp_type == "experimental" & df$dpi > 0] = "after" # change the class label for experimental samples after infection
df$dol = as.numeric(df$dol) # change dol to a numeric
pre_axis_cutoff_df = df

# Plot HM days
ggplot(pre_axis_cutoff_df, aes(x=Axis.1, y=Axis.2, color=HM_days)) +
  geom_point(size=2) +
  # scale_color_manual(values=c("black", "red")) +
  scale_color_gradient(low="orange", high="blue") +
  theme(
    axis.text = element_text(size=25),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  xlab(paste0("PCoA1 (", pc1, ")")) +
  ylab(paste0("PCoA2 (", pc2, ")"))


# Plot DOL
ggplot(pre_axis_cutoff_df, aes(x=Axis.1, y=Axis.2, color=dol)) +
  geom_point(size=2) +
  # scale_color_manual(values=c("black", "red")) +
  scale_color_gradient(low="orange", high="blue") +
  theme(
    axis.text = element_text(size=25),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  xlab(paste0("PCoA1 (", pc1, ")")) +
  ylab(paste0("PCoA2 (", pc2, ")"))


# Plot Formula days
ggplot(pre_axis_cutoff_df, aes(x=Axis.1, y=Axis.2, color=Form_days)) +
  geom_point(size=2) +
  # scale_color_manual(values=c("black", "red")) +
  scale_color_gradient(low="orange", high="blue") +
  theme(
    axis.text = element_text(size=25),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  xlab(paste0("PCoA1 (", pc1, ")")) +
  ylab(paste0("PCoA2 (", pc2, ")"))

