setwd("/Users/drewschwartz/Box Sync/Neonate BSI/Manuscript/Tables")
metadata<-read.csv("MasterMetadata_FINAL_GitHub.csv", sep=",") # make sure you save the metadata as a csv
uniq_subjects = metadata[!(duplicated(metadata$Subject)),]
controls = uniq_subjects[uniq_subjects$exp_type == "control",]
exps = uniq_subjects[uniq_subjects$exp_type == "experimental",]
write.csv(uniq_subjects, "unique.csv")
Allcontrols = metadata[which(metadata$exp_type %in% c("control")),]
Allexpts = metadata[which(metadata$exp_type %in% c("experimental")),]
#did the above to determine that there are 210 experimental stools and 252 control stools

#factor = "abx_score_complete" # change this variable for whichever variable you're investigating

#wilcox.test(controls[, c(factor)], exps[, c(factor)])
#median(controls[, c(factor)])
#median(exps[, c(factor)])
#IQR(controls[, c(factor)])
#IQR(exps[, c(factor)])

#median(controls[, c(factor)])
#[1] 129
# median(exps[, c(factor)])
#[1] 142
# IQR(controls[, c(factor)])
#[1] 80
# IQR(exps[, c(factor)])
#[1] 108
#p = 0.3031

factor1 = "gest_age" # change this variable for whichever variable you're investigating

wilcox.test(controls[, c(factor1)], exps[, c(factor1)])
median(controls[, c(factor1)])
median(exps[, c(factor1)])
IQR(controls[, c(factor1)])
IQR(exps[, c(factor1)])
quantile(controls[, c(factor1)])
quantile(exps[, c(factor1)])
#median(controls[, c(factor1)])
#[1] 25
#median(exps[, c(factor1)])
#[1] 25
#IQR(controls[, c(factor1)])
#1
#IQR(exps[, c(factor1)])
#[1] 3
# quantile(controls[, c(factor1)])
# 0%  25%  50%  75% 100% 
# 23   25   25   26   29 
# > quantile(exps[, c(factor1)])
# 0%  25%  50%  75% 100% 
# 23   24   25   27   29 
#p = 0.6473

factor2 = "birthweight" # change this variable for whichever variable you're investigating

wilcox.test(controls[, c(factor2)], exps[, c(factor2)])
median(controls[, c(factor2)])
median(exps[, c(factor2)])
IQR(controls[, c(factor2)])
IQR(exps[, c(factor2)])
quantile(controls[, c(factor2)])
quantile(exps[, c(factor2)])
# Controls 810 IQR204; Experimental 740 IQR212.5
#p = 0.2223
# quantile(controls[, c(factor2)])
# 0%  25%  50%  75% 100% 
# 400  676  810  880 1100 
# > quantile(exps[, c(factor2)])
# 0%    25%    50%    75%   100% 
# 400.0  620.0  740.0  832.5 1041.0 

#Controls: 6/37 vaginal = 16.2%; 31/37 caesarean = 83.7%
#Cases: 3/19 vaginal = 15.8%; 16/19 caesarean = 84%
#Controls: 12/37 male = 32.4%; 25/37 female = 67.8%
#Cases: 6/19 male = 31.5%; 13/19 female = 68.4%

factor3 = "culture_bug1_DOL" # change this variable for whichever variable you're investigating

wilcox.test(controls[, c(factor3)], exps[, c(factor3)])
median(controls[, c(factor3)])
median(exps[, c(factor3)])
IQR(controls[, c(factor3)])
IQR(exps[, c(factor3)])
quantile(controls[, c(factor3)])
quantile(exps[, c(factor3)])
#Bacteremia DOL Controls: 32 IQR23; Cases: 32 IQR24.5; p = 0.9103
# quantile(controls[, c(factor3)])
# 0%  25%  50%  75% 100% 
# 13   19   32   46   81 
# > quantile(exps[, c(factor3)])
# 0%  25%  50%  75% 100% 
# 13.0 19.5 32.0 43.5 81.0 

factor4 = "ABx_days_beforebacteremia"
wilcox.test(controls[, c(factor4)], exps[, c(factor4)])
median(controls[, c(factor4)])
median(exps[, c(factor4)])
IQR(controls[, c(factor4)])
IQR(exps[, c(factor4)])
quantile(controls[, c(factor4)])
quantile(exps[, c(factor4)])
#Controls 21 IQR 19
#Cases 16 IQR 18
#p = 0.7032
# quantile(controls[, c(factor4)])
# 0%  25%  50%  75% 100% 
# 0    9   21   28   50 
# > quantile(exps[, c(factor4)])
# 0%  25%  50%  75% 100% 
# 0    7   16   25   43 

factor5 = "ABx_score_beforebacteremia" # change this variable for whichever variable you're investigating

wilcox.test(controls[, c(factor5)], exps[, c(factor5)])
median(controls[, c(factor5)])
median(exps[, c(factor5)])
IQR(controls[, c(factor5)])
IQR(exps[, c(factor5)])
quantile(controls[, c(factor5)])
quantile(exps[, c(factor5)])
#Controls 57 IQR 60
#Case 53 IQR 86
#p = 0.5796
# quantile(controls[, c(factor5)])
# 0%  25%  50%  75% 100% 
# 0   25   57   85  122 
# > quantile(exps[, c(factor5)])
# 0%  25%  50%  75% 100% 
# 0   26   53  112  255 

Allcontrolsbefore<-subset(Allcontrols, DPI1<=0)
Allcontrolsafter<-subset(Allcontrols, DPI1>0)
Allexptsbefore<-subset(Allexpts, DPI1<=0)
Allexptsafter<-subset(Allexpts, DPI1>0)


