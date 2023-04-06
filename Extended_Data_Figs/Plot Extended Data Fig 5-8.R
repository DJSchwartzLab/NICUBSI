library(ggplot2)
library(RColorBrewer)
library(reshape2)


setwd("/path/to/github/repo/NICUBSI/tables")

# GENERATE METAPHLAN LINE/DOT PLOTS

species = read.table("211210_wide_metaphlanStandard_species_100k.txt", header=T, sep=" ")
species$dol = sapply(strsplit(as.character(species$Sample), "_"), function(s) s[2])
species$id = sapply(strsplit(as.character(species$Sample), "_"), function(s) s[1])

exp_ids = c("130.01", "137.01", "145.01", "1117.01","1057.01","1002.01", "396.01", "362.01", "340.01","297.01",
            "273.01", "254.01", "225.01", "205.01", "96.01", "87.01","86.01", "151.01", "175.01")
bugs = c("s__Enterococcus_faecalis", "s__Serratia_marcescens", "s__Escherichia_coli", "s__Enterobacter_cloacae_complex", 
         "s__Klebsiella_pneumoniae", "s__Enterococcus_faecalis","s__Staphylococcus_aureus","s__Staphylococcus_aureus",
         "s__Klebsiella_pneumoniae", "s__Escherichia_coli", "s__Serratia_marcescens", "s__Staphylococcus_aureus",
         "s__Staphylococcus_aureus", "s__Serratia_marcescens", "s__Staphylococcus_aureus", "s__Streptococcus_agalactiae",
         "s__Staphylococcus_aureus", "s__Streptococcus_agalactiae", "s__Streptococcus_agalactiae")
infx_dols = c(45, 13, 20, 34, 35, 42, 18, 74, 48, 32, 19, 21, 14, 14, 27, 50, 26, 41, 81)

exp_data = data.frame(ids=exp_ids, bug=bugs, infx_dol=infx_dols)
opath = "/path/to/output_dir"

dir.create(opath)

for(exp_id in as.character(exp_data$ids)){
  subject_data = species[species$id == exp_id,]
  subject_data$dol = gsub("a|b|c", "", subject_data$dol)
  bug = exp_data[match(exp_id, exp_data$ids), c("bug")]
  infx_dol = exp_data[match(exp_id, exp_data$ids), c("infx_dol")]
  sub_data = subject_data[, c("dol", as.character(bug))]
  colnames(sub_data) = c("dol","bug")
  ggplot(sub_data, aes(x=as.numeric(dol), y=bug)) + 
    ylim(limits=c(0, 100)) +
    geom_vline(xintercept=infx_dol, color="red", linetype="dashed", size=1) +
    geom_line(size=1) +
    geom_point(size=3) +
    theme(legend.position = "none",
          axis.text = element_text(size=25),
          axis.title = element_text(size=20),
          panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.grid.major = element_line(color="grey95", size=0.7), 
          panel.grid.minor = element_line(color="grey95", size=0.7),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    xlab("DOL") +
    ggtitle(exp_id) +
    ylab(paste(substr(as.character(bug), 4, nchar(as.character(bug))), "Rel. Abundance"))
  ggsave(filename=paste0(exp_id, "_0_100_v3.pdf"), path=opath, units="in", width=12, height=7)
}

