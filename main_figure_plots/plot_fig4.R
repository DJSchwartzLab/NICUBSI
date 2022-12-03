library(reshape2)
library(vegan)
library(ecodist)
library(ape)
library(plyr)
library(dplyr)
library(colortools)
library(data.table)
library(ggplot2)
library(lmerTest)
library(lme4)
library(lsmeans)

# set the current working dir to the path of the GitHub repository
setwd("/path/to/github/repo/NICUBSI/tables")

# PLOT 3A
metadata = read.csv("MasterMetadata_FINAL_GitHub.csv")

species = read.table("211210_wide_metaphlanStandard_species_100k.txt", header=T, sep=" ")
species$id = sapply(strsplit(as.character(species$Sample), "_"), function(s) s[1])
species$dol = sapply(strsplit(as.character(species$Sample), "_"), function(s) s[2])

# collate the relevant metadata for the metaphlan calls
species$exp_type = metadata$exp_type[match(species$Sample, metadata$Sequence.name)]
species$exp_id = metadata$exp_id[match(species$Sample, metadata$Sequence.name)]
species$DPI = metadata$DPI1[match(species$Sample, metadata$Sequence.name)]
species = species[species$exp_type %in% c("control", "experimental"), ]

exp_ids = c("130.01", "137.01", "145.01", "1117.01","1057.01","1002.01", "396.01", "362.01", "340.01","297.01",
            "273.01", "254.01", "225.01", "205.01", "96.01", "87.01","86.01", "151.01", "175.01")

# hard code metaphlan names for pathogen species
bugs = c("s__Enterococcus_faecalis", "s__Serratia_marcescens", "s__Escherichia_coli", "s__Enterobacter_cloacae_complex", 
         "s__Klebsiella_pneumoniae", "s__Enterococcus_faecalis","s__Staphylococcus_aureus","s__Staphylococcus_aureus",
         "s__Klebsiella_pneumoniae", "s__Escherichia_coli", "s__Serratia_marcescens", "s__Staphylococcus_aureus",
         "s__Staphylococcus_aureus", "s__Serratia_marcescens", "s__Staphylococcus_aureus", "s__Streptococcus_agalactiae",
         "s__Staphylococcus_aureus", "s__Streptococcus_agalactiae", "s__Streptococcus_agalactiae"
)
# SANITY CHECK! 
paste(bugs, metadata$culture_bug1[match(exp_ids, metadata$Subject)])

exp_data = data.frame(ids=exp_ids, bug=bugs)
species$bug_select = sapply(1:dim(species)[1], function(i){ # THIS IS FOR 2A
  if(species[i, c("exp_type")] == "control"){
    "none"
  }else{
    bug = exp_data$bug[match(species[i, c("id")], exp_data$ids)][1]
    bug_abundance = species[i, c(as.character(bug))] # pathogen we care about
    # if(is.null(bug_abundance)) bug_abundance = 0
    bug_distrib = as.vector(species[i, grepl("s__", colnames(species))]) # all species counts in sample
    1 - (sum(bug_distrib > bug_abundance) / length(bug_distrib[bug_distrib > 0]))
  }
})

case_subjects = species[species$exp_type == "experimental",]

paste(case_subjects$DPI, cut(case_subjects$DPI, breaks=c(-15, -10, -5, 0, 5, 10, 15), right=F))
case_subjects$dpi_range = as.character(cut(case_subjects$DPI, breaks=c(-15, -10, -5, 0, 5, 10, 15), right=F))
case_subjects = case_subjects[!is.na(case_subjects$dpi_range),]
case_subjects$dpi_range
first_5 = case_subjects[case_subjects$dpi_range == "[-15,-10)",]
second_5 = case_subjects[case_subjects$dpi_range == "[-10,-5)",]
third_5 = case_subjects[case_subjects$dpi_range == "[-5,0)",]
fourth_5 = case_subjects[case_subjects$dpi_range == "[0,5)",]
fifth_5 = case_subjects[case_subjects$dpi_range == "[5,10)",]
sixth_5 = case_subjects[case_subjects$dpi_range == "[10,15)",]


first_5_df = data.frame(ids=unique(first_5$id), 
                        med_rank=sapply(unique(first_5$id), function(s) median(as.numeric(unlist(first_5$bug_select[first_5$id == s])))),
                        time_range=c("-15 to -10"))
second_5_df = data.frame(ids=unique(second_5$id), 
                         med_rank=sapply(unique(second_5$id), function(s) median(as.numeric(unlist(second_5$bug_select[second_5$id == s])))),
                         time_range=c("-10 to -5"))
third_5_df = data.frame(ids=unique(third_5$id), 
                        med_rank=sapply(unique(third_5$id), function(s) median(as.numeric(unlist(third_5$bug_select[third_5$id == s])))),
                        time_range=c("-5 to 0"))
fourth_5_df = data.frame(ids=unique(fourth_5$id), 
                         med_rank=sapply(unique(fourth_5$id), function(s) median(as.numeric(unlist(fourth_5$bug_select[fourth_5$id == s])))),
                         time_range=c("0 to 5"))
fifth_5_df = data.frame(ids=unique(fifth_5$id), 
                        med_rank=sapply(unique(fifth_5$id), function(s) median(as.numeric(unlist(fifth_5$bug_select[fifth_5$id == s])))),
                        time_range=c("5 to 10"))
sixth_5_df = data.frame(ids=unique(sixth_5$id), 
                        med_rank=sapply(unique(sixth_5$id), function(s) median(as.numeric(unlist(sixth_5$bug_select[sixth_5$id == s])))),
                        time_range=c("10 to 15"))

bin_ranks = function(rank){
  rank = as.numeric(as.character(rank))
  if(is.na(rank)){
    rank = 0
  }
  if(rank >= .75){
    "100-75"
  }else if(rank >= .50){
    "74-50"
  }else if(rank >= .25){
    "49-25"
  }else{
    "24-0"
  }
}

first_5_df$rank_bin = sapply(first_5_df$med_rank, bin_ranks)
second_5_df$rank_bin = sapply(second_5_df$med_rank, bin_ranks)
third_5_df$rank_bin = sapply(third_5_df$med_rank, bin_ranks)
fourth_5_df$rank_bin = sapply(fourth_5_df$med_rank, bin_ranks)
fifth_5_df$rank_bin = sapply(fifth_5_df$med_rank, bin_ranks)
sixth_5_df$rank_bin = sapply(sixth_5_df$med_rank, bin_ranks)

bin_labels = c("100-75", "74-50", "49-25", "24-0")
first_5_freqs=data.frame(bin=bin_labels, range="-15 to -10", freq=sapply(bin_labels, function(bin) dim(first_5_df[first_5_df$rank_bin == bin,])[1] / dim(first_5_df)[1]))
second_5_freqs=data.frame(bin=bin_labels, range="-10 to -5",freq=sapply(bin_labels, function(bin) dim(second_5_df[second_5_df$rank_bin == bin,])[1] / dim(second_5_df)[1]))
third_5_freqs=data.frame(bin=bin_labels, range="-5 to 0",freq=sapply(bin_labels, function(bin) dim(third_5_df[third_5_df$rank_bin == bin,])[1] / dim(third_5_df)[1]))
fourth_5_freqs=data.frame(bin=bin_labels, range="0 to 5",freq=sapply(bin_labels, function(bin) dim(fourth_5_df[fourth_5_df$rank_bin == bin,])[1] / dim(fourth_5_df)[1]))
fifth_5_freqs=data.frame(bin=bin_labels, range="5 to 10",freq=sapply(bin_labels, function(bin) dim(fifth_5_df[fifth_5_df$rank_bin == bin,])[1] / dim(fifth_5_df)[1]))
sixth_5_freqs=data.frame(bin=bin_labels, range="10 to 15",freq=sapply(bin_labels, function(bin) dim(sixth_5_df[sixth_5_df$rank_bin == bin,])[1] / dim(sixth_5_df)[1]))
full_df = rbind(first_5_freqs, second_5_freqs, third_5_freqs, fourth_5_freqs, fifth_5_freqs, sixth_5_freqs)
full_df$range = factor(full_df$range, levels=c("-15 to -10", "-10 to -5", "-5 to 0", "0 to 5", "5 to 10", "10 to 15"))
full_df$bin = factor(full_df$bin, levels=c("100-75", "74-50","49-25", "24-0"))

# 3a
ggplot(full_df, aes(x=range, y=freq)) +
  geom_bar(aes(fill=bin), color="black", width=0.9, stat="identity", position="stack") +
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.0))) +
  scale_x_discrete(expand = expansion(mult = c(0.0, 0.0))) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "YlGnBu")) +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

# run chi-squared 
cast_full_df = data.table::dcast(data=full_df, formula=bin ~ range)
rownames(cast_full_df) = cast_full_df$bin
cast_full_df = cast_full_df[, 2:ncol(cast_full_df)]
chisq.test(cast_full_df)

# plot and compute stats for 2b and 2c

instrain_classes = read.delim("220424_genome_summary_removedunder100K.tsv")

instrain_classes$subject = sapply(strsplit(as.character(instrain_classes$sample), "_"), function(s) paste(s[[1]][1], sep = "_")) # parse sample names
instrain_classes$strain = sapply(strsplit(as.character(instrain_classes$sample), "_"), function(s) paste(s[grep("DJS", s):length(s)], collapse="_"))
instrain_classes$dol = sapply(strsplit(as.character(instrain_classes$sample), "_"), function(s) paste(s[[2]][1], sep = "_"))
# remove trailing characters in dols
instrain_classes$dol_parse = sapply(instrain_classes$dol, function(s){
  if(grepl("a|b", s)){
    substr(s, 1, nchar(s)-1)
  }else{
    s
  }
})

instrain_classes$dol_parse = as.numeric(as.character(instrain_classes$dol_parse))
instrain_classes$dol_infx = metadata$culture_bug1_DOL[match(instrain_classes$subject, metadata$Subject)]
instrain_classes$DPI = instrain_classes$dol_parse - as.numeric(as.character(instrain_classes$dol_infx))
# create binning functions
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

instrain_classes$seq_name = paste(instrain_classes$subject, instrain_classes$dol, sep="_")
instrain_classes$read_count = metadata$Post.processed_reads[match(instrain_classes$seq_name, metadata$Sequence.name)]
instrain_classes$read_count = as.numeric(as.character(instrain_classes$read_count))
instrain_classes$sample[is.na(instrain_classes$read_count)]

instrain_classes$dpi_bin = sapply(instrain_classes$DPI, bin)

subj_max = function(s){
  max(instrain_classes$cov_normalized[instrain_classes$subject == s])
}
subj_min = function(s){
  min(instrain_classes$cov_normalized[instrain_classes$subject == s])
}
instrain_classes = instrain_classes[instrain_classes$read_count > 1e5 & instrain_classes$dpi_bin !="NA",]
instrain_classes$group = paste(instrain_classes$subject, instrain_classes$dpi_bin)
instrain_classes$cov_normalized = instrain_classes$coverage_median / (instrain_classes$read_count) # normalize isolate coverage by read counts
instrain_classes$cov_normalized_0to1 = (instrain_classes$cov_normalized - sapply(instrain_classes$subject, subj_min)) / sapply(instrain_classes$subject, subj_max)

# linear mixed effect model for stats
cov_normalized_stats = lmer(cov_normalized_0to1~DPI+(1|Subject), data=instrain_classes[instrain_classes$DPI >= 0,])
anova(cov_normalized_stats)
post_hoc_corrected_stats = lsmeans(cov_normalized_stats, cov_normalized_0to1~DPI, adjust="BH")

# plot 3B
ggplot(instrain_classes, aes(x=DPI, y=cov_normalized_0to1)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_smooth(formula = y ~ x, se=T, color="red", size=1.5, level=0.95, fill="lightcoral") +
  geom_point(alpha=0.2, color="grey50") +
  scale_x_continuous(limits=c(-12, NA), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits=c(NA, 1.02), expand = expansion(mult = c(0, 0))) +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

########################################################################################################################################################
# 3D
########################################################################################################################################################

genome_summary = read.delim("220424_genome_summary_removedunder100K.tsv")

genome_summary$id = sapply(strsplit(as.character(genome_summary$sample), "_"), function(s) s[1] )


genome_summary$bug = tolower(sapply(strsplit(as.character(genome_summary$sample), "_"), function(s) paste(s[grep("DJS.*", s) + 2], s[grep("DJS.*", s) + 3], sep="_") )) # parse sample names for the pathogen otu
genome_summary$bug[grepl("mrsa", genome_summary$bug)] = "staphylococcus_aureus" # fix pathogen names
genome_summary$bug[grepl("klebsiella_pneumonia", genome_summary$bug)] = "klebsiella_pneumoniae" # fix pathogen names

genome_summary$log_sns = -log10(genome_summary$parse_sns/(genome_summary$length*genome_summary$breadth))
genome_summary$log_breadth = -log10(1-genome_summary$breadth)
genome_summary$log_sns[is.infinite(genome_summary$log_sns)] = 7 # for plotting purposes so that perfect cases appear on plot range
genome_summary$log_breadth[is.infinite(genome_summary$log_breadth)] = 7 # for plotting purposes so that perfect cases appear on plot range
genome_summary$sns_rate = genome_summary$parse_sns/(genome_summary$length*genome_summary$breadth)

ggplot(genome_summary[genome_summary$breadth > 0.1,], aes(x=log_sns, y=log_breadth)) + # plot scatter
  geom_point(aes(color=bug, size=log_sns)) +
  scale_size_continuous(range=c(0.25, 2.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits=c(1, NA)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits=c(1, NA)) +
  theme(
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

ggplot(genome_summary[genome_summary$breadth > 0.1 ,], aes(x=log_sns)) + # plot first adjoined histogram in 2d
  geom_histogram(aes(y = (..count..)/sum(..count..)), fill="lightskyblue3", color="grey1", bins=10) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits=c(1, 7), breaks=1:7, labels = as.character(seq(1,7))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

ggplot(genome_summary[genome_summary$breadth > 0.1 ,], aes(x=log_breadth)) + # plot second adjoined histogram in 2d
  geom_histogram(aes(y = (..count..)/sum(..count..)), fill="lightskyblue3", color="grey1", bins=10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), position="right") +
  scale_x_reverse(expand = expansion(mult = c(0.02, 0)), limits=c(7,1), breaks=seq(7,1), labels = as.character(seq(7,1))) +
  theme(
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

########################################################################################################################################################
# 3E
########################################################################################################################################################

genome_summary$id = sapply(strsplit(as.character(genome_summary$sample), "_"), function(s) s[1] )
genome_summary_fill = genome_summary[!is.na(genome_summary$parse_sns),]

genome_summary_fill$dol = sapply(strsplit(as.character(genome_summary_fill$sample), "_"), function(s) s[2] ) # parse out dol from sample name
genome_summary_fill$dol = sub("a|b|c", "", genome_summary_fill$dol) # remove trailinig characters
genome_summary_fill$day_of_infection = metadata$culture_bug1_DOL[match(genome_summary_fill$id, metadata$Subject)] # map the day of life infection to the subject
genome_summary_fill$dpi = as.numeric(as.character(genome_summary_fill$dol)) - genome_summary_fill$day_of_infection

# construct plottable df
sns_all = data.frame(id=unique(genome_summary_fill$id), 
                     sns_counts=as.character(sapply(unique(genome_summary_fill$id), function(id) min(genome_summary_fill$parse_sns[genome_summary_fill$id == id & genome_summary_fill$breadth > 0.99]))))

sns_all$sns_counts = as.numeric(as.character(sns_all$sns_counts))

# sns counts before infx
sns_all$sns_counts_before_infx = as.character(sapply(unique(genome_summary_fill$id), function(id){ 
  min(genome_summary_fill$parse_sns[genome_summary_fill$id == id & 
                                      genome_summary_fill$dpi <= 0 & 
                                      genome_summary_fill$breadth > 0.99])
  
}))
sns_all_keep = unique(sns_all$id[sns_all$sns_counts < 20])
sns_df = data.frame(sns_num=as.character(c(sns_all$sns_counts_before_infx, sns_all$sns_counts)), type=c(rep("before infx", times=dim(sns_all)[1]), rep("all", times=dim(sns_all)[1])) )
sns_df$sns_num = as.numeric(as.character(sns_df$sns_num))

ggplot(sns_df, aes(x=sns_num)) + # 0 - 10 SNS
  geom_bar(aes(fill=type,), stat="count", position=position_dodge()) +
  scale_x_continuous(limits=c(NA, 10.5), breaks=0:10, expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks=0:10, expand = expansion(mult = c(0, 0.05))) +
  theme(
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major.y = element_line(color="grey95", size=0.7), 
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "black", size=0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values=c("coral3", "cornflowerblue")) +
  labs(fill = "date range")

ggplot(sns_df, aes(x=sns_num)) + # 0 - 20 SNS
  geom_bar(aes(fill=type,), stat="count", position=position_dodge()) +
  scale_x_continuous(limits=c(NA, 20.5), breaks=0:20, expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks=0:10, expand = expansion(mult = c(0, 0.05))) +
  theme(
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major.y = element_line(color="grey95", size=0.7), 
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "black", size=0.5),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values=c("coral3", "cornflowerblue")) +
  labs(fill = "date range")


########################################################################################################################################################
# 3C
########################################################################################################################################################

# filter plotted subjects to only those who display sns < 20 at 1+ timepoints
instrain_classes = instrain_classes[instrain_classes$subject %in% sns_all_keep,]

df_list = lapply(unique(as.character(instrain_classes$group)), function(g) {
  df = instrain_classes[instrain_classes$group == g,]
  df$cov_normalized = mean(as.numeric(as.character(unlist(df$cov_normalized))))
  df$reads_mean = mean(as.numeric(as.character(unlist(df$read_count))))
  df[1 ,]
})

lumped_coverage_df = as.data.frame(rbindlist(df_list))
lumped_coverage_df$dpi_bin = factor(lumped_coverage_df$dpi_bin, levels=c("-15 to -10", 
                                                                         "-10 to -5", 
                                                                         "-5 to 0", 
                                                                         "0 to 5", 
                                                                         "5 to 10", 
                                                                         "10 to 15"))

# compute wilcoxon p values between DPI groups
combo_bin = t(as.data.frame(combn(levels(lumped_coverage_df$dpi_bin), m=2)))
colnames(combo_bin) = c("DPI1", "DPI2")
combo_bin = as.data.frame(combo_bin)
combo_bin$wilcox_p = sapply(1:nrow(combo_bin), function(i){
  dpi_1_group = combo_bin$DPI1[i]
  dpi_2_group = combo_bin$DPI2[i]
  wilcox.res = wilcox.test(
    lumped_coverage_df$cov_normalized[lumped_coverage_df$dpi_bin == dpi_1_group], lumped_coverage_df$cov_normalized[lumped_coverage_df$dpi_bin == dpi_2_group]
  )
  wilcox.res$p.value
})

# BH correction
combo_bin$corrected_wilcox_p = pmin(1, combo_bin$wilcox_p * nrow(combo_bin))
write.table(combo_bin, file="3C_pairwise_wilcox.tsv", sep="\t", row.names=F)

# Plot 3C
ggplot(lumped_coverage_df, aes(x=dpi_bin, y=cov_normalized)) +
  geom_boxplot(fill="peachpuff2") +
  geom_jitter(width=0.05, alpha=0.4) +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=20),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

