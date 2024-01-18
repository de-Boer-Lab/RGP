#SUPP FIG 13
z_score_1 = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/CAGE_analysis/z_filt_data_1.csv', sep = ",", header = TRUE)
## Wilcox test
wilcox.test(z_score_1$genomic, z_score_1$di.local)

z_score_3 = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/CAGE_analysis/z_filt_data_3.csv', sep = ",", header = TRUE)

## Wilcox test
wilcox.test(z_score_3$genomic, z_score_3$di.local)


z_score_10 = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/CAGE_analysis/z_filt_data_10.csv', sep = ",", header = TRUE)

## Wilcox test
wilcox.test(z_score_10$genomic, z_score_10$di.local)

z_score_50 = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/CAGE_analysis/z_filt_data_50.csv', sep = ",", header = TRUE)
## Wilcox test
wilcox.test(z_score_50$genomic, z_score_50$di.local)

z_score_100 = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/CAGE_analysis/z_filt_data_100.csv', sep = ",", header = TRUE)

## Wilcox test
wilcox.test(z_score_100$genomic, z_score_100$di.local)

#FIG 4 B
meso_endo_ecto = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/meso_endo_ecto/meso_endo_ecto.csv', sep = ",", header = TRUE)
cor.test(meso_endo_ecto$Mesoderm.naive, meso_endo_ecto$Mesoderm.evolved, alternative = c("greater"))
cor.test(meso_endo_ecto$Mesoderm.naive, meso_endo_ecto$Endoderm.naive)
cor.test(meso_endo_ecto$Mesoderm.naive, meso_endo_ecto$Ectoderm.naive)

cor.test(meso_endo_ecto$Mesoderm.naive, meso_endo_ecto$Mesoderm.evolved)
cor.test(meso_endo_ecto$Endoderm.naive, meso_endo_ecto$Endoderm.evolved)
cor.test(meso_endo_ecto$Endoderm.naive, meso_endo_ecto$Ectoderm.naive)

cor.test(meso_endo_ecto$Ectoderm.evolved, meso_endo_ecto$Mesoderm.evolved)
cor.test(meso_endo_ecto$Ectoderm.evolved, meso_endo_ecto$Endoderm.evolved)
cor.test(meso_endo_ecto$Ectoderm.naive, meso_endo_ecto$Ectoderm.evolved)

#SUPP FIG 13
dnase_inpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/Dnase_inpeak.csv', sep = ",", header = TRUE)
dnase_notinpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/Dnase_notinpeak.csv', sep = ",", header = TRUE)
wilcox.test(dnase_inpeak$Dnase, dnase_notinpeak$Dnase)

H3K4me3_inpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K4me3_inpeak.csv', sep = ",", header = TRUE)
H3K4me3_notinpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K4me3_notinpeak.csv', sep = ",", header = TRUE)
wilcox.test(H3K4me3_inpeak$H3K4me3, H3K4me3_notinpeak$H3K4me3)


H3K27me3_inpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K27me3_inpeak.csv', sep = ",", header = TRUE)
H3K27me3_notinpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K27me3_notinpeak.csv', sep = ",", header = TRUE)
wilcox.test(H3K27me3_inpeak$H3K27me3, H3K27me3_notinpeak$H3K27me3)


H3K4me1_inpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K4me1_inpeak.csv', sep = ",", header = TRUE)
H3K4me1_notinpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K4me1_notinpeak.csv', sep = ",", header = TRUE)
wilcox.test(H3K4me1_inpeak$H3K4me1, H3K4me1_notinpeak$H3K4me1)

H3K27ac_inpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K27ac_inpeak.csv', sep = ",", header = TRUE)
H3K27ac_notinpeak = read.table('/Users/ishika/Desktop/PhD_Project/Enformer_Analysis/Paper_Review_Analysis/Analysis_output_06152023/H1_DNase_rep1_07052023/Evolution_plots/H3K27ac_notinpeak.csv', sep = ",", header = TRUE)
wilcox.test(H3K27ac_inpeak$H3K27ac, H3K27ac_notinpeak$H3K27ac)
