library(optparse)
library(dplyr)
library(ggplot2)
library(plyr)
library(data.table)


# set the current working dir to the path of the GitHub repository
setwd("/path/to/github/repo/NICUBSI/tables")

# PLOT 2a
metadata = read.csv("MasterMetadata_FINAL_GitHub.csv")
species = read.table("211210_wide_metaphlanStandard_species_100k.txt", header=T, sep=" ")

rownames(species) = species$Sample
species = species[, grepl("s__", colnames(species))]
diversity = data.frame(Sample=rownames(species), shannon=as.numeric(vegan::diversity(species, "shannon")), simpson=as.numeric(vegan::diversity(species, "simpson")))

diversity$DPI = master_data$DPI1[match(diversity$Sample, master_data$Sequence.name)]
diversity$exp = master_data$exp_type[match(diversity$Sample, master_data$Sequence.name)]
diversity$bug = master_data$culture_bug1[match(diversity$Sample, master_data$Sequence.name)]

diversity$id = sapply(strsplit(as.character(diversity$Sample), "_"), function(c) c[[1]][1])
diversity$exp = factor(diversity$exp, levels=c("control", "experimental"))
diversity = diversity[diversity$exp %in% c("control", "experimental") & !is.na(diversity$DPI) ,]

bin = function(d){
  if(abs(d) > 15){
    "NA"
  }else if(d <= -10){
    "-15 to -10"
  }else if(d <= -5){
    "-10 to -5"
  }else if(d <= 0){
    "-5 to 0"
  }else if(d <= 5){
    "0 to 5"
  }else if(d <= 10){
    "5 to 10"
  }else{
    "10 to 15"
  }
}

diversity$dpi_class = sapply(diversity$DPI, function(d) bin(d))
diversity$group = paste(diversity$dpi_class, diversity$id)
df_list = lapply(unique(as.character(diversity$group)), function(g) {
  df = diversity[diversity$group == g,]
  df$shannon = mean(df$shannon)
  df$simpson = mean(df$simpson)
  df[1 ,]
})

lumped_div_df = rbindlist(df_list)
lumped_div_df = lumped_div_df[!(lumped_div_df$dpi_class == "NA")]
lumped_div_df$dpi_class = factor(lumped_div_df$dpi_class, levels=c("-15 to -10", "-10 to -5", "-5 to 0",
                                                                   "0 to 5", "5 to 10", "10 to 15"))
# Plot 2A
ggplot(lumped_div_df, aes(x=dpi_class, y=shannon, fill=exp)) +
  geom_boxplot(position=position_dodge(0.7), width=0.5, outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(0.1), alpha=0.4) +
  theme(
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  scale_fill_manual(values=c("gray26", "red3"))


# Plot 2C
metadata = read.csv("MasterMetadata_FINAL_GitHub.csv")
species = read.table("211210_wide_metaphlanStandard_species_100k.txt", header=T, sep=" ")
species$id = sapply(strsplit(as.character(species$Sample), "_"), function(s) s[1])
species$dol = sapply(strsplit(as.character(species$Sample), "_"), function(s) s[2])

# collate the relevant metadata for the metaphlan calls
species$exp_type = metadata$exp_type[match(species$Sample, metadata$Sequence.name)]
species$exp_id = metadata$exp_id[match(species$Sample, metadata$Sequence.name)]

species$exp_id[grepl("68.01", species$Sample)]
species$DPI = metadata$DPI1[match(species$Sample, metadata$Sequence.name)]
species = species[species$exp_type %in% c("control", "experimental"), ]

exp_ids = c("130.01", "137.01", "145.01", "1117.01","1057.01","1002.01", "396.01", "362.01", "340.01","297.01",
            "273.01", "254.01", "225.01", "205.01", "96.01", "87.01","86.01", "151.01", "175.01")

# hard code metaphlan names for pathogen species
bugs = c("s__Enterococcus_faecalis", "s__Serratia_marcescens", "s__Escherichia_coli", "s__Enterobacter_cloacae", 
         "s__Klebsiella_pneumoniae", "s__Enterococcus_faecalis","s__Staphylococcus_aureus","s__Staphylococcus_aureus",
         "s__Klebsiella_pneumoniae", "s__Escherichia_coli", "s__Serratia_marcescens", "s__Staphylococcus_aureus",
         "s__Staphylococcus_aureus", "s__Serratia_marcescens", "s__Staphylococcus_aureus", "s__Streptococcus_agalactiae",
         "s__Staphylococcus_aureus", "s__Streptococcus_agalactiae", "s__Streptococcus_agalactiae"
)
# SANITY CHECK! AND I AM SANE, FOR THE MOMENT
paste(bugs, metadata$culture_bug1[match(exp_ids, metadata$Subject)])

species$bug_select = sapply(1:dim(species)[1], function(i){ # retrieving bug relative abundance
  if(species[i, c("exp_type")] == "control"){
    bug = exp_data$bug[match(species[i, c("exp_id")], exp_data$ids)][1]
    species[i, c(as.character(bug))]
  }else{
    bug = exp_data$bug[match(species[i, c("id")], exp_data$ids)][1]
    species[i, c(as.character(bug))]
  }
})


bin = function(d){
  if(abs(d) > 14){
    "NA"
  }else{
    "-14 to 14"
  }
}

species = species[!is.na(species$DPI),]
species = species[species$DPI <= 0,]
species$dpi_class = sapply(as.numeric(as.character(species$DPI)), function(d) bin(d))

species$group = paste(species$dpi_class, species$id)

df_list = lapply(unique(as.character(species$group)), function(g) {
  df = species[species$group == g,]
  df$bug_select = mean(as.numeric(as.character(unlist(df$bug_select))))
  df[1 ,]
})

lumped_div_df = as.data.frame(rbindlist(df_list))

lumped_div_df = lumped_div_df[!(as.character(lumped_div_df$dpi_class) == "NA"),]
lumped_div_df$dpi_class = factor(lumped_div_df$dpi_class, levels=c( "-14 to 14"))

sapply(levels(lumped_div_df$dpi_class), function(c) wilcox.test(lumped_div_df$bug_select[lumped_div_df$dpi_class == c & lumped_div_df$exp_type == "control"], 
                                                                lumped_div_df$bug_select[lumped_div_df$dpi_class == c & lumped_div_df$exp_type == "experimental"], alternative="l"))

wilcox.test(lumped_div_df$bug_select[lumped_div_df$exp_type == "control"], lumped_div_df$bug_select[lumped_div_df$exp_type == "experimental"])

control_quantile = paste(round(quantile(lumped_div_df$bug_select[!is.na(lumped_div_df$dpi_class) & lumped_div_df$exp_type == "control"], probs=c(0.25, 0.5, 0.75), na.rm=T), 2), collapse=",")
exp_quantile = paste(round(quantile(lumped_div_df$bug_select[!is.na(lumped_div_df$dpi_class) & lumped_div_df$exp_type == "experimental"], probs=c(0.25, 0.5, 0.75), na.rm=T), 2), collapse=",")

ggplot(lumped_div_df[!is.na(lumped_div_df$dpi_class),], aes(x=exp_type, y=bug_select, fill=exp_type)) +
  geom_boxplot(position=position_dodge(0.7), outlier.shape = NA) +
  geom_jitter(width=0.05, alpha=0.5) +
  scale_fill_manual(values=c("black", "red")) +
  ggtitle(paste("control", control_quantile, "experimental", exp_quantile)) +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


# PLOT 2D

# load metadata
metadata = read.csv("MasterMetadata_FINAL_GitHub.csv")

# load the species relative abundance

# sample specific factors: Abx_days, ABx_score_sample, DOL, NPO_days, HM_days, Form_days,	DM_days, DPI1
# subject specific factors: exp_type, gest_age, birthweight, sex, bug1_Family

species_df = read.delim("211210_wide_metaphlanStandard_species_100k.txt", sep=" ")

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
species_mat = read.delim("211210_wide_metaphlanStandard_species_100k.txt", sep=" ")
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
df = df[, !grepl("Axis", colnames(df))]
rownames(df) = df$sample_name
rownames(df) = 1:nrow(df)

permanove_distance_matrix = vegdist(sapply(species, as.numeric),  method = "bray", na.rm=T)

subject_types = df$subject
sample_types = 1:length(df$subject)
non_duplicated_subjects = !duplicated(subject_types)
# sample specific factors: Abx_days, ABx_score_sample, DOL, NPO_days, HM_days, Form_days,	DM_days, DPI1
# subject specific factors: exp_type, gest_age, birthweight, sex, bug1_Family

# sample specific df
abx_df = data.frame(abx=df$abx_score)
abx_days_df = data.frame(abx_sample_days=df$abx_days)
dol_df = data.frame(DOL=df$dol)
NPO_days_df = data.frame(npo_days=df$NPO_days)
HM_days_df = data.frame(hm_days=df$HM_days)
Form_days_df = data.frame(form_days=df$Form_days)
dpi_df = data.frame(dpi=df$dpi)

# subject specific df
exp_type_df = data.frame(exp_type=df$exp_type[non_duplicated_subjects], row.names=unique(subject_types))
ga_df = data.frame(ga=df$ga[non_duplicated_subjects], row.names=unique(subject_types))
bw_df = data.frame(bw=df$bw[non_duplicated_subjects], row.names=unique(subject_types))
sex_df = data.frame(sex=df$sex[non_duplicated_subjects], row.names=unique(subject_types))
bug_family_df = data.frame(bug_family=df$Bug1_Family[non_duplicated_subjects], row.names=unique(subject_types))
subject_df = data.frame(subject=df$subject[non_duplicated_subjects], row.names=unique(subject_types))

# permanova call - subject specific
abx_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  sample_data = abx_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)

abx_days_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  sample_data = abx_days_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)

dol_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  sample_data = dol_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)

npo_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  sample_data = NPO_days_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)

hm_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  sample_data = HM_days_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
) # plot the PCOA

form_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  sample_data = Form_days_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
) # plot the PCOA

dpi_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  sample_data = dpi_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
) # plot the PCOA

# permanova call - subject specific
exp_type_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  subject_data = exp_type_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)
ga_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  subject_data = ga_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)
bw_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  subject_data = bw_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)
sex_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  subject_data = sex_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)
bug_family_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  subject_data = bug_family_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)
subject_results = PERMANOVA_repeat_measures(
  permanove_distance_matrix, # create a bray-curtis distance matrix
  subject_types, # ordered vector with subject id matching samples
  subject_data = subject_df, # ordered table with metadata factors we want to compute significance for
  permutations = 999,
  ncores = 1
)

# p_raw = c(0.805, 0.793, 0.001, 0.801, 0.001, 0.019, 0.474, 0.163, 0.789, 0.868, 0.092, 0.001, 0.002)
# p_adj = p.adjust(p_raw, method="BH")

ggplot(pre_axis_cutoff_df, aes(x=Axis.1, y=Axis.2, color=Bug1_Family)) +
  geom_point(size=2) +
  # scale_color_manual(values=c("black", "red")) +
  # scale_color_gradient(low="orange", high="blue") +
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


