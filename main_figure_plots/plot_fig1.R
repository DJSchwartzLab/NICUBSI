input_data_familyDGMG<-read.delim("Metagenomics/MollyDrewSRAMetaphlan3/220301_wide_metaphlanStandard_family.txt", header = TRUE)
metadataDGMGtest<-read.delim("Metadata/220508_gaspgibsonSRANICU_metadata.csv", sep=",", header=TRUE)
row.names(metadataDGMGtest)<-metadataDGMGtest$Sample

row.names(input_data_familyDGMG)<-input_data_familyDGMG$Sample
input_data_familyDGMG$Sample=NULL

fit_data = Maaslin2(
  input_data=input_data_familyDGMG, input_metadata=metadataDGMG, output= "Maaslin/DrewGMollyGMeta3/last10d_familymicrobiomeampgentvanc",
  fixed_effects = c("Ampicillin_10d","Gentamicin_10d","Vancomycin_10d","AmpxGent","VancxGent","DOL"), random_effects = c("Subject_id"), transform = "AST",reference = c('AmpxGent,neither','VancxGent,neither'))

