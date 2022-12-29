####PLOT EXTENDED DATA FIGURE 1B###

#Antibiotic score before bacteremia by DPI bin
ggplot(metadata, aes(x=dpi_bin, y=ABx_score_sample, fill=exp_type)) +
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
  scale_fill_manual(values=c("gray26", "red3"))+scale_x_discrete(limits=c("-15 to -10", "-10 to -5", "-5 to 0", "0 to 5", "5 to 10", "10 to 15"))

#15d before and after
abxscore15<-subset(metadata, DPI1<=15)
abxscore15<-subset(abxscore15, DPI1>=-15)
abxscore15$dpi_binexp<-paste(abxscore15$exp_type,abxscore15$dpi_bin)
fitlogABXscoreBIN<- lmer(ABx_score_sample~dpi_binexp+(1|Subject), data=abxscore15)
anova(fitlogABXscoreBIN)
# dpi_binexp   11  59892  5444.7  44.651
tukfitlogABXscoreBIN<-lsmeans(fitlogABXscoreBIN, pairwise~dpi_binexp, adjust="tukey")
tukfitlogABXscoreBIN
# contrast                                             estimate    SE    df t.ratio p.value
#  (control -10 to -5) - (control -15 to -10)               4.03  2.70 333.8   1.492  0.9421
#  (control -10 to -5) - (control -5 to 0)                 -7.06  2.58 333.4  -2.730  0.2163
#  (control -10 to -5) - control 0 to 5                   -10.56  2.54 333.8  -4.152  0.0024
#  (control -10 to -5) - control 10 to 15                 -21.90  2.75 333.4  -7.951  <.0001
#  (control -10 to -5) - control 5 to 10                  -17.00  2.78 334.4  -6.123  <.0001
#  (control -10 to -5) - (experimental -10 to -5)         -23.63 13.23  60.9  -1.786  0.8190
#  (control -10 to -5) - (experimental -15 to -10)        -14.34 13.29  62.0  -1.079  0.9946
#  (control -10 to -5) - (experimental -5 to 0)           -24.83 13.02  57.1  -1.908  0.7496
#  (control -10 to -5) - experimental 0 to 5              -39.59 13.04  57.5  -3.036  0.1247
#  (control -10 to -5) - experimental 10 to 15            -65.44 13.23  60.8  -4.947  0.0004
#  (control -10 to -5) - experimental 5 to 10             -60.08 13.07  58.1  -4.595  0.0013
#  (control -15 to -10) - (control -5 to 0)               -11.08  2.57 333.5  -4.319  0.0012
#  (control -15 to -10) - control 0 to 5                  -14.59  2.58 334.3  -5.651  <.0001
#  (control -15 to -10) - control 10 to 15                -25.92  2.85 334.2  -9.104  <.0001
#  (control -15 to -10) - control 5 to 10                 -21.03  2.82 334.9  -7.452  <.0001
#  (control -15 to -10) - (experimental -10 to -5)        -27.65 13.23  60.8  -2.091  0.6310
#  (control -15 to -10) - (experimental -15 to -10)       -18.37 13.29  62.0  -1.382  0.9628
#  (control -15 to -10) - (experimental -5 to 0)          -28.86 13.02  57.1  -2.217  0.5447
#  (control -15 to -10) - experimental 0 to 5             -43.62 13.04  57.5  -3.345  0.0586
#  (control -15 to -10) - experimental 10 to 15           -69.47 13.23  60.8  -5.252  0.0001
#  (control -15 to -10) - experimental 5 to 10            -64.11 13.07  58.1  -4.903  0.0005
#  (control -5 to 0) - control 0 to 5                      -3.50  2.46 334.0  -1.427  0.9576
#  (control -5 to 0) - control 10 to 15                   -14.84  2.68 333.6  -5.545  <.0001
#  (control -5 to 0) - control 5 to 10                     -9.94  2.70 334.6  -3.678  0.0143
#  (control -5 to 0) - (experimental -10 to -5)           -16.57 13.21  60.5  -1.255  0.9817
#  (control -5 to 0) - (experimental -15 to -10)           -7.28 13.27  61.6  -0.549  1.0000
#  (control -5 to 0) - (experimental -5 to 0)             -17.77 12.99  56.7  -1.368  0.9650
#  (control -5 to 0) - experimental 0 to 5                -32.53 13.02  57.1  -2.499  0.3616
#  (control -5 to 0) - experimental 10 to 15              -58.38 13.21  60.5  -4.421  0.0023
#  (control -5 to 0) - experimental 5 to 10               -53.02 13.05  57.7  -4.062  0.0075
#  control 0 to 5 - control 10 to 15                      -11.34  2.62 333.9  -4.327  0.0012
#  control 0 to 5 - control 5 to 10                        -6.44  2.64 334.8  -2.436  0.3852
#  control 0 to 5 - (experimental -10 to -5)              -13.07 13.19  60.2  -0.991  0.9974
#  control 0 to 5 - (experimental -15 to -10)              -3.78 13.26  61.4  -0.285  1.0000
#  control 0 to 5 - (experimental -5 to 0)                -14.27 12.98  56.4  -1.100  0.9935
#  control 0 to 5 - experimental 0 to 5                   -29.03 13.00  56.9  -2.233  0.5342
#  control 0 to 5 - experimental 10 to 15                 -54.88 13.19  60.2  -4.161  0.0053
#  control 0 to 5 - experimental 5 to 10                  -49.52 13.04  57.5  -3.798  0.0167
#  control 10 to 15 - control 5 to 10                       4.90  2.85 334.4   1.716  0.8598
#  control 10 to 15 - (experimental -10 to -5)             -1.73 13.24  61.1  -0.131  1.0000
#  control 10 to 15 - (experimental -15 to -10)             7.56 13.31  62.3   0.568  1.0000
#  control 10 to 15 - (experimental -5 to 0)               -2.93 13.03  57.3  -0.225  1.0000
#  control 10 to 15 - experimental 0 to 5                 -17.69 13.06  57.8  -1.355  0.9673
#  control 10 to 15 - experimental 10 to 15               -43.54 13.24  61.1  -3.288  0.0665
#  control 10 to 15 - experimental 5 to 10                -38.18 13.09  58.4  -2.917  0.1619
#  control 5 to 10 - (experimental -10 to -5)              -6.63 13.23  61.0  -0.501  1.0000
#  control 5 to 10 - (experimental -15 to -10)              2.66 13.30  62.1   0.200  1.0000
#  control 5 to 10 - (experimental -5 to 0)                -7.83 13.02  57.2  -0.601  1.0000
#  control 5 to 10 - experimental 0 to 5                  -22.59 13.05  57.6  -1.732  0.8460
#  control 5 to 10 - experimental 10 to 15                -48.44 13.23  61.0  -3.661  0.0242
#  control 5 to 10 - experimental 5 to 10                 -43.08 13.08  58.2  -3.294  0.0666
#  (experimental -10 to -5) - (experimental -15 to -10)     9.29  4.19 333.4   2.218  0.5378
#  (experimental -10 to -5) - (experimental -5 to 0)       -1.20  3.22 333.3  -0.374  1.0000
#  (experimental -10 to -5) - experimental 0 to 5         -15.96  3.37 333.4  -4.741  0.0002
#  (experimental -10 to -5) - experimental 10 to 15       -41.81  3.95 333.3 -10.598  <.0001
#  (experimental -10 to -5) - experimental 5 to 10        -36.46  3.40 333.1 -10.712  <.0001
#  (experimental -15 to -10) - (experimental -5 to 0)     -10.49  3.44 333.4  -3.049  0.0994
#  (experimental -15 to -10) - experimental 0 to 5        -25.25  3.64 333.7  -6.939  <.0001
#  (experimental -15 to -10) - experimental 10 to 15      -51.10  4.28 333.7 -11.950  <.0001
#  (experimental -15 to -10) - experimental 5 to 10       -45.74  3.63 333.3 -12.589  <.0001
#  (experimental -5 to 0) - experimental 0 to 5           -14.76  2.46 333.6  -6.009  <.0001
#  (experimental -5 to 0) - experimental 10 to 15         -40.61  3.28 333.7 -12.382  <.0001
#  (experimental -5 to 0) - experimental 5 to 10          -35.25  2.54 333.4 -13.899  <.0001
#  experimental 0 to 5 - experimental 10 to 15            -25.85  3.31 333.5  -7.809  <.0001
#  experimental 0 to 5 - experimental 5 to 10             -20.49  2.66 333.5  -7.702  <.0001
#  experimental 10 to 15 - experimental 5 to 10             5.36  3.43 333.4   1.563  0.9211

###Extended Data Figure 1C Plotting###

#Shortbred Plotting! This first section is to determine the median #of ARGs in the 14d prior to BSI for cases vs. controls
sample_df = data.frame(samples=unique(abx_genes$sample), count=sapply(unique(abx_genes$sample), function(s) dim(abx_genes[abx_genes$sample == s,])[1]))

sample_df$dpi = metadata$DPI1[match(sample_df$samples, metadata$Sequence.name)]
sample_df$subject = metadata$Subject[match(sample_df$samples, metadata$Sequence.name)]
sample_df_dpi_restricted = sample_df[sample_df$dpi <= 0 & sample_df$dpi >= -14,]

subj_collapsed_df = aggregate(sample_df_dpi_restricted$count, list(sample_df_dpi_restricted$subject), median)
colnames(subj_collapsed_df) = c("subject", "gene_count")
subj_collapsed_df$exp_type = metadata$exp_type[match(subj_collapsed_df$subject, metadata$Subject)]

w_test = wilcox.test(subj_collapsed_df$gene_count[subj_collapsed_df$exp_type == "control"], subj_collapsed_df$gene_count[subj_collapsed_df$exp_type == "experimental"])
control_iqr = paste(round(quantile(subj_collapsed_df$gene_count[subj_collapsed_df$exp_type == "control"], c(0.25, 0.5, 0.75), na.rm=T), 1), collapse=",")
exp_iqr = paste(round(quantile(subj_collapsed_df$gene_count[subj_collapsed_df$exp_type == "experimental"], c(0.25, 0.5, 0.75), na.rm=T), 1), collapse=",")
ggplot(subj_collapsed_df, aes(y=gene_count, x=exp_type, fill=exp_type)) +
  geom_boxplot() +
  geom_jitter(width=0.05, alpha=0.5) +scale_fill_manual(values=c("gray26", "red3"))+
  ggtitle(paste0("dpi<0, p=", round(w_test$p.value, 3), control_iqr, exp_iqr)) +
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=20),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_line(color="grey95", size=0.7), 
    panel.grid.minor = element_line(color="grey95", size=0.7),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
