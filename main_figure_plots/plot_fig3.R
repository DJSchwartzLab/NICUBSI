#Figure 3A
genus_taxa <- read.delim("Taxa_contribution_for_R.txt")
metadata<-read.delim("MasterMetadata_FINAL_GitHub.csv")
#set color palette for bar plots
palette_all <- grDevices::colors()
palette_no_grey <- palette_all[grep("gr(a|e)y",
                                    grDevices::colors(),
                                    invert=T)]
n_colors <- 100
set.seed(267)
my_palette <- sample(palette_no_grey,n_colors)

metadata<-subset(metadata,Post.processed_reads>100000)
#Any difference in enterotype assignment by bug1_family for all data
tblgenus=table(metadata$Bug1_Family, metadata$Genus_enterotype_assignment)
chisq.test(tblgenus, correct = FALSE)
#X-squared = 106.84, df = 12, p-value < 2.2e-16
tblgenus
#                     m1  m2  m3  m4
#  Control            107  55  44  38
#  Enterobacteriaceae  60  26   5   7
#  Enterococcaceae     17   2   0  13
#  Staphylococcaceae   40   1  10   0
#  Streptococcaceae     7   1  15   1

tblpairGenus<-pairwiseNominalIndependence(tblgenus,
                                          fisher = FALSE,
                                          gtest  = FALSE,
                                          chisq  = TRUE,
                                          method = "fdr")
view(tblpairGenus)

genus_cont %>% 
  ggplot(aes(Model,Abundance,fill=Taxa))+
  geom_bar(stat="identity",position="stack")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = my_palette)#+


#Figures 3B-F

#Combining our dataframe and Molly/Drew dataframes (run before changing sample to rownames)
input_data_species <- read.delim("211210_wide_metaphlanStandard_species_100k.txt", header = TRUE) #The abundance table file
input_data_genus <- read.delim("211210_wide_metaphlanStandard_genus.txt")
# The metadata table file
input_metadata <- read.delim("MasterMetadata_FINAL_GitHub.csv", header = TRUE, sep = ",") # The metadata table file
input_metadata[input_metadata=="N/A"] <-NA
input_metadata<-subset(input_metadata,Post.processed_reads>100000)
metadataDGMG<-read.delim("221102_DrewMollyMeta.csv", sep=",", header=TRUE)

combinedspecies<-rbind.fill(input_data_species,input_data_speciesDGMG)
combinedspecies[is.na(combinedspecies)]<-0
row.names(combinedspecies)<-combinedspecies$Sample
combinedspecies$Sample=NULL

combinedgenus<-rbind.fill(input_data_genus, input_data_genusDGMG)
combinedgenus[is.na(combinedgenus)]<-0
row.names(combinedgenus)<-combinedgenus$Sample
combinedgenus$Sample=NULL

#Now combine metadata
input_metadata$Sample<-input_metadata$Sequence.name
combinedmetadata<-rbind.fill(input_metadata, metadataDGMG)
row.names(combinedmetadata)<-combinedmetadata$Sample

#Now, combined Maaslin for pre-bacteremia versus all others controlled for DOL!!!
#Enterococcus
preEnterococcuscombined<-subset(combinedmetadata, Enterococcus_bacteremia=="yes")
preEnterococcuscombined<-subset(preEnterococcuscombined, DPI1<=0)
noEnterococcuscombined<-subset(combinedmetadata, Enterococcus_bacteremia=="no")

testEnterococcuscombined<-rbind(preEnterococcuscombined,noEnterococcuscombined)
testEnterococcuscombinedspecies<-merge(testEnterococcuscombined, combinedspecies, by='row.names')
# write.csv(testEnterococcuscombinedspecies, "Metagenomics/MollyDrewNICUBSIcombined/221215_NICUBSItestEnterococcuscombined.csv")
fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=testEnterococcuscombined, output= "Maaslin/NICUBSIbacteremia/221214preEnterococcusCombined",
  fixed_effects = c("Enterococcus_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")

#S. Aureus
preStaphcombined<-subset(combinedmetadata, StaphylococcusAureus_bacteremia=="yes")
preStaphcombined<-subset(preStaphcombined, DPI1<=0)
noStaphcombined<-subset(combinedmetadata, StaphylococcusAureus_bacteremia=="no")

testStaphcombined<-rbind(preStaphcombined,noStaphcombined)
testStaphcombinedspecies<-merge(testStaphcombined, combinedspecies,by='row.names')
# write.csv(testStaphcombinedspecies, "Metagenomics/MollyDrewNICUBSIcombined/221215_NICUBSItestStaphaureuscombined.csv")
fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=testStaphcombined, output= "Maaslin/NICUBSIbacteremia/221214preStaphaureusCombined",
  fixed_effects = c("StaphylococcusAureus_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")

testStaphcombined<-rbind(preStaphcombined,noStaphcombined)
fit_data = Maaslin2(
  input_data=combinedgenus, input_metadata=testStaphcombined, output= "Maaslin/NICUBSIbacteremia/221214preStaphaureusCombinedgenus",
  fixed_effects = c("Staphylococcus_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")

#GBS
preGBScombined<-subset(combinedmetadata, Streptococcus_bacteremia=="yes")
preGBScombined<-subset(preGBScombined, DPI1<=0)
noGBScombined<-subset(combinedmetadata, Streptococcus_bacteremia=="no")

testGBScombined<-rbind(preGBScombined,noGBScombined)
testGBScombinedspecies<-merge(testGBScombined, combinedspecies,by='row.names')
# write.csv(testGBScombinedspecies, "Metagenomics/MollyDrewNICUBSIcombined/221215_NICUBSItestGBScombined.csv")
fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=testGBScombined, output= "Maaslin/NICUBSIbacteremia/221214preGBSCombinedminprev0.05",
  fixed_effects = c("Streptococcus_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST", min_prevalence = 0.05)
fit_data = Maaslin2(
  input_data=combinedgenus, input_metadata=testGBScombined, output= "Maaslin/NICUBSIbacteremia/preGBSCombinedgenus",
  fixed_effects = c("Streptococcus_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")

#All Enterobacteriaceae
preEnterobacteriaceaecombined<-subset(combinedmetadata, Enterobacteriaceae_bacteremia=="yes")
preEnterobacteriaceaecombined<-subset(preEnterobacteriaceaecombined, DPI1<=0)
noEnterobacteriaceacombined<-subset(combinedmetadata, Enterobacteriaceae_bacteremia=="no")
testEnterobacteriaceaecombined<-rbind(preEnterobacteriaceaecombined,noEnterobacteriaceacombined)

fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=testEnterobacteriaceaecombined, output= "Maaslin/NICUBSIbacteremia/221214preEnterobacteriaceaeCombined",
  fixed_effects = c("Enterobacteriaceae_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")
fit_data = Maaslin2(
  input_data=combinedgenus, input_metadata=testEnterobacteriaceaecombined, output= "Maaslin/NICUBSIbacteremia/221214preEnterobacteriaceaeCombinedgenus",
  fixed_effects = c("Enterobacteriaceae_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")
fit_data = Maaslin2(
  input_data=combinedfamily, input_metadata=testEnterobacteriaceaecombined, output= "Maaslin/NICUBSIbacteremia/preEnterobacteriaceaeCombinedfamily",
  fixed_effects = c("Enterobacteriaceae_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")

#Serratia
preSerratiacombined<-subset(combinedmetadata, Serratia_bacteremia=="yes")
preSerratiacombined<-subset(preSerratiacombined, DPI1<=0)
noSerratiacombined<-subset(combinedmetadata, Serratia_bacteremia=="no")
testSerratiacombined<-rbind(preSerratiacombined,noSerratiacombined)
testSerratiaombinedspecies<-merge(testSerratiacombined, combinedspecies,by='row.names')
# write.csv(testSerratiaombinedspecies, "Metagenomics/MollyDrewNICUBSIcombined/221215_NICUBSItestSerratiacombined.csv")
#Since only 50 kids had Serratia at all in the gut microbiome, I had to reduce the min_prevalance feature from 0.1 default to 0.5
fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=testSerratiacombined, output= "Maaslin/NICUBSIbacteremia/221214preSerratiaCombined",
  fixed_effects = c("Serratia_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST", min_prevalence = 0.05)

fit_data = Maaslin2(
  input_data=combinedgenus, input_metadata=testSerratiacombined, output= "Maaslin/NICUBSIbacteremia/preSerratiaCombinedgenus",
  fixed_effects = c("Serratia_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST", min_prevalence = 0.05)

#Klebsiella
preKlebsiellacombined<-subset(combinedmetadata, Klebpneumo_bacteremia=="yes")
preKlebsiellacombined<-subset(preKlebsiellacombined, DPI1<=0)
noKlebsiellacombined<-subset(combinedmetadata, Klebpneumo_bacteremia=="no")
testKlebsiellacombined<-rbind(preKlebsiellacombined,noKlebsiellacombined)
testKlebsiellacombinedgenus<-merge(testKlebsiellacombined, combinedgenus,by='row.names')
# write.csv(testKlebsiellacombinedgenus, "Metagenomics/MollyDrewNICUBSIcombined/221215_NICUBSItestKlebsiellaacombined.csv")
fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=testKlebsiellacombined, output= "Maaslin/NICUBSIbacteremia/221214preKlebCombined",
  fixed_effects = c("Klebpneumo_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")
fit_data = Maaslin2(
  input_data=combinedgenus, input_metadata=testKlebsiellacombined, output= "Maaslin/NICUBSIbacteremia/221214preKlebCombinedgenus",
  fixed_effects = c("Klebpneumo_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")

#E. coli
preEcolicombined<-subset(combinedmetadata, Ecoli_bacteremia=="yes")
preEcolicombined<-subset(preEcolicombined, DPI1<=0)
noEcolicombined<-subset(combinedmetadata, Ecoli_bacteremia=="no")
testEcolicombined<-rbind(preEcolicombined,noEcolicombined)
testEcolicombinedspecies<-merge(testEcolicombined, combinedspecies,by='row.names')
# write.csv(testEcolicombinedspecies, "Metagenomics/MollyDrewNICUBSIcombined/221215_NICUBSItestecolicombined.csv")
fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=testEcolicombined, output= "Maaslin/NICUBSIbacteremia/221214preEcoliCombined",
  fixed_effects = c("Ecoli_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST")

#Testing whether all by all changes findings.
fit_data = Maaslin2(
  input_data=combinedspecies, input_metadata=combinedmetadata, output= "Maaslin/NICUBSIbacteremia/ALLCOMBINEDSPECIES",
  fixed_effects = c("Enterococcus_bacteremia","Ecoli_bacteremia","Klebpneumo_bacteremia","Serratia_bacteremia","Streptococcus_bacteremia","Staphylococcus_bacteremia","DOL"), random_effects = c("Subject"), transform = "AST", min_prevalence = 0.05)
