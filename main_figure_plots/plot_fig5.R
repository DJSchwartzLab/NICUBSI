library(reshape2)
library(vegan)
library(ecodist)
library(ape)
library(plyr)
library(dplyr)
library(colortools)
library(data.table)
library(ggplot2)

# set the current working dir to the path of the GitHub repository
setwd("/path/to/github/repo/NICUBSI/tables")

# load case mapping
case_mapping = read.csv("metagenomesmappedCaseisolates_COMBINED.csv")

# plot 5A

case_mapping$negative_log_breadth = -log10(1 - case_mapping$breadth)
case_mapping$negative_log_breadth[is.infinite(case_mapping$negative_log_breadth)] = 7
case_mapping$negative_log_sns_rate = -log10(case_mapping$SNS_count / case_mapping$length)
case_mapping$negative_log_sns_rate[is.infinite(case_mapping$negative_log_sns_rate)] = 7

write.csv(case_mapping, file="metagenomesmappedCaseisolates_COMBINED_fix.csv", quote=F, row.names=F)

case_mapping$case_to_case = case_mapping$metagenome_subject == case_mapping$Isolate_subject

ggplot(case_mapping, aes(x=negative_log_sns_rate, y=negative_log_breadth,)) +
  geom_point(aes(color=case_to_case), size=2.5) +
  scale_color_manual(values=c("black", "red")) +
  scale_y_continuous(limits=c(0.5, 7), expand = expansion(mult = c(0.01, 0.05)), breaks = c(0.5, 2:7), labels = c(0.5, 2:7)) +
  scale_x_continuous(limits=c(1.5, 7), expand = expansion(mult = c(0.01, 0.05)), breaks = c(1.5, 2:7), labels = c(1.5, 2:7)) + 
  theme(
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# plot 5B

case_mapping$breadth = as.numeric(as.character(case_mapping$breadth))
case_mapping$SNS_count = as.numeric(as.character(case_mapping$SNS_count))
case_mapping$length = as.numeric(as.character(case_mapping$length))

case_mapping$breadth_sns_stat_raw = case_mapping$SNS_count  / (case_mapping$breadth_minCov * case_mapping$length)
# if above is 0 -> SNS_count is 0
# if above is NA -> some NA in the raw inputs (leave alone)
case_mapping$breadth_sns_stat = -log10(case_mapping$breadth_sns_stat_raw)
case_mapping$breadth_sns_stat[is.infinite(case_mapping$breadth_sns_stat) & (case_mapping$breadth_sns_stat > 0)] = 7
case_mapping$breadth_sns_stat[is.infinite(case_mapping$breadth_sns_stat) & (case_mapping$breadth_sns_stat < 0)] = 0

case_mapping$species = sapply(strsplit(case_mapping$sample, "_"), function(sample_split){
  vec_len = length(sample_split)
  if("MRSA" %in% sample_split){
    start = vec_len
  }else{
    start = vec_len - 1
  }
  species = sample_split[start: vec_len]
  species = paste0(species, collapse="_")
  tolower(species)
})

species_name_fix_df = data.frame(original=unique(case_mapping$species), fix=c(
  "enterococcus_faecalis", "klebsiella_pneumoniae", "enterobacter_cloacae", "klebsiella_pneumoniae", "staphylococcus_aureus", "serratia_marcescens", "escherichia_coli", "serratia_marcescens", "group_b_steptococcus", "staphylococcus_aureus")
)
# Enterococcaceae, Staphylococcaceae, Streptococcaceae, Enterobacteriales
species_name_fix_df$family = c("Enterococcaceae", "Enterobacteriales", "Enterobacteriales", "Enterobacteriales", "Staphylococcaceae", "Enterobacteriales", "Enterobacteriales", "Enterobacteriales", "Streptococcaceae", "Staphylococcaceae")

case_mapping$species_fix = species_name_fix_df$fix[match(case_mapping$species, species_name_fix_df$original)]
case_mapping$family_fix = species_name_fix_df$family[match(case_mapping$species, species_name_fix_df$original)]

match_case = match(case_mapping$species, species_name_fix_df$original)

case_mapping$mg_subject_isolate_key = paste(case_mapping$metagenome_subject, case_mapping$Isolate_subject, sep="_") 
case_mapping_order = case_mapping[order(case_mapping$breadth_sns_stat, decreasing=T),]
case_mapping_unique = case_mapping[!duplicated(case_mapping$mg_subject_isolate_key),]

# PLOT 5B
ggplot(case_mapping_unique[case_mapping_unique$breadth > 0.5 & case_mapping_unique$case_to_case == F,], aes(x=family_fix, y=breadth_sns_stat,)) +
  geom_jitter(height=0, width=0.2, aes(color=breadth_sns_stat > 5)) +
  scale_color_manual(values=c("black", "#2b95ff")) +
  geom_vline(xintercept=15.5, linetype="dashed", color="grey50") +
  scale_y_continuous(limits=c(2, 7), expand = expansion(mult = c(0.01, 0.05)), breaks = c(0.5, 2:7), labels = c(0.5, 2:7)) +
  # scale_x_continuous(limits=c(0, 7), expand = expansion(mult = c(0.01, 0.05))) + 
  theme(
    axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.3),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

