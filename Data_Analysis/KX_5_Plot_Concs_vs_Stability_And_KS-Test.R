# Create plot for concentration vs stability

conc_file = "Data/Metabolite_Concentration_Sets.txt"
perc_file_Base = "Data/KX_Results_1_met_set_vs_percent_steady_Base.tab"
perc_file_FSBPase = "Data/KX_Results_1_met_set_vs_percent_steady_FSBPase.tab"
perc_file_TKT = "Data/KX_Results_1_met_set_vs_percent_steady_TKT.tab"
perc_file_BOTH = "Data/KX_Results_1_met_set_vs_percent_steady_BOTH.tab"

library(reshape2)
library(tidyverse)
library(foreach)
library(doMC)
library(ggridges)
library(scales)
library(optparse)
library(ggrepel)
library(ggplot2)

# Load data
conc_data = read.table(conc_file, header=T, sep="\t")
conc_data = conc_data[1:5000,]
perc_data_Base = read.table(perc_file_Base, header=F, sep="\t")
perc_data_FSBPase = read.table(perc_file_FSBPase, header=F, sep="\t")
perc_data_TKT = read.table(perc_file_TKT, header=F, sep="\t")
perc_data_BOTH = read.table(perc_file_BOTH, header=F, sep="\t")

colnames(perc_data_Base) = c("set", "stable")
colnames(perc_data_FSBPase) = c("set", "stable")
colnames(perc_data_TKT) = c("set", "stable")
colnames(perc_data_BOTH) = c("set", "stable")


# Transform the data
conc_data = log10(conc_data)

# Split data into lower and upper stable steady state categories
perc_data_Base = perc_data_Base[order(perc_data_Base$stable, decreasing=T),]
perc_data_FSBPase = perc_data_FSBPase[order(perc_data_FSBPase$stable, decreasing=T),]
perc_data_TKT = perc_data_TKT[order(perc_data_TKT$stable, decreasing=T),]
perc_data_BOTH = perc_data_BOTH[order(perc_data_BOTH$stable, decreasing=T),]

# Divide metabolite concentration sets into high and low % stable steady states


# By Deciles
deciles_Base = quantile(perc_data_Base$stable, prob = seq(0, 1, length = 11), type = 5)
d10_Base = deciles_Base[10]
d1_Base = deciles_Base[2]

perc_hi_Base = subset(perc_data_Base, stable > d10_Base)
perc_lo_Base = subset(perc_data_Base, stable <= d1_Base)

# Subset concentrations to those in high and low % stable steady states groups
conc_hi_Base = conc_data[perc_hi_Base$set,]
conc_lo_Base = conc_data[perc_lo_Base$set,]




deciles_FSBPase = quantile(perc_data_FSBPase$stable, prob = seq(0, 1, length = 11), type = 5)
d10_FSBPase = deciles_FSBPase[10]
d1_FSBPase = deciles_FSBPase[2]

perc_hi_FSBPase = subset(perc_data_FSBPase, stable > d10_FSBPase)
perc_lo_FSBPase = subset(perc_data_FSBPase, stable <= d1_FSBPase)

# Subset concentrations to those in high and low % stable steady states groups
conc_hi_FSBPase = conc_data[perc_hi_FSBPase$set,]
conc_lo_FSBPase = conc_data[perc_lo_FSBPase$set,]




deciles_TKT = quantile(perc_data_TKT$stable, prob = seq(0, 1, length = 11), type = 5)
d10_TKT = deciles_TKT[10]
d1_TKT = deciles_TKT[2]

perc_hi_TKT = subset(perc_data_TKT, stable > d10_TKT)
perc_lo_TKT = subset(perc_data_TKT, stable <= d1_TKT)

# Subset concentrations to those in high and low % stable steady states groups
conc_hi_TKT = conc_data[perc_hi_TKT$set,]
conc_lo_TKT = conc_data[perc_lo_TKT$set,]



deciles_BOTH = quantile(perc_data_BOTH$stable, prob = seq(0, 1, length = 11), type = 5)
d10_BOTH = deciles_BOTH[10]
d1_BOTH = deciles_BOTH[2]

perc_hi_BOTH = subset(perc_data_BOTH, stable > d10_BOTH)
perc_lo_BOTH = subset(perc_data_BOTH, stable <= d1_BOTH)

# Subset concentrations to those in high and low % stable steady states groups
conc_hi_BOTH = conc_data[perc_hi_BOTH$set,]
conc_lo_BOTH = conc_data[perc_lo_BOTH$set,]





# Plot it
library(ggplot2)

## Base
c_lo_long_Base = reshape2::melt(conc_lo_Base)
colnames(c_lo_long_Base) = c("metabolite", "concentration")
c_lo_long_Base$stability = rep("unstable", nrow(c_lo_long_Base))

c_hi_long_Base = reshape2::melt(conc_hi_Base)
colnames(c_hi_long_Base) = c("metabolite", "concentration")
c_hi_long_Base$stability = rep("stable", nrow(c_hi_long_Base))

c_long_Base = rbind(c_hi_long_Base, c_lo_long_Base)
c_long_Base$variant = rep("4_Base", nrow(c_long_Base))

## FSBPase
c_lo_long_FSBPase = reshape2::melt(conc_lo_FSBPase)
colnames(c_lo_long_FSBPase) = c("metabolite", "concentration")
c_lo_long_FSBPase$stability = rep("unstable", nrow(c_lo_long_FSBPase))

c_hi_long_FSBPase= reshape2::melt(conc_hi_FSBPase)
colnames(c_hi_long_FSBPase) = c("metabolite", "concentration")
c_hi_long_FSBPase$stability = rep("stable", nrow(c_hi_long_FSBPase))

c_long_FSBPase = rbind(c_hi_long_FSBPase, c_lo_long_FSBPase)
c_long_FSBPase$variant = rep("3_FSBPase", nrow(c_long_FSBPase))

## TKT
c_lo_long_TKT = reshape2::melt(conc_lo_TKT)
colnames(c_lo_long_TKT) = c("metabolite", "concentration")
c_lo_long_TKT$stability = rep("unstable", nrow(c_lo_long_TKT))

c_hi_long_TKT = reshape2::melt(conc_hi_TKT)
colnames(c_hi_long_TKT) = c("metabolite", "concentration")
c_hi_long_TKT$stability = rep("stable", nrow(c_hi_long_TKT))

c_long_TKT = rbind(c_hi_long_TKT, c_lo_long_TKT)
c_long_TKT$variant = rep("2_TKT", nrow(c_long_TKT))

## BOTH
c_lo_long_BOTH = reshape2::melt(conc_lo_BOTH)
colnames(c_lo_long_BOTH) = c("metabolite", "concentration")
c_lo_long_BOTH$stability = rep("unstable", nrow(c_lo_long_BOTH))

c_hi_long_BOTH = reshape2::melt(conc_hi_BOTH)
colnames(c_hi_long_BOTH) = c("metabolite", "concentration")
c_hi_long_BOTH$stability = rep("stable", nrow(c_hi_long_BOTH))

c_long_BOTH = rbind(c_hi_long_BOTH, c_lo_long_BOTH)
c_long_BOTH$variant = rep("1_BOTH", nrow(c_long_BOTH))



c_long_All = rbind(c_long_Base,c_long_FSBPase,c_long_TKT,c_long_BOTH)
MetOfInterest = c("Xu5P","RuBP","GAP","F6P","E4P","S7P","R5P","SBP","Ru5P","AMP")
c_long_Sel = c_long_All[c_long_All$metabolite %in% MetOfInterest, ]





c_long_All_MOD = c_long_All
Mets_to_Remove = c("CO2_cax","CO2_cyt","O2")#,"FUMPool","COAPool","PPool")
for(i in 1:length(Mets_to_Remove)){
  c_long_All_MOD = subset(c_long_All_MOD, metabolite!=Mets_to_Remove[i])
}



ggplot(c_long_All_MOD, aes(x = concentration, y = variant, fill = stability)) +
  #stat_density_ridges(scale = 0.85, alpha = 0.5,quantile_lines = TRUE, quantiles = 0.5) +
  geom_density_ridges(scale = 0.85, alpha = 0.5) +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme_bw() +
  scale_fill_manual(values=c("#8073ac","#e08214")) +
  theme(strip.background = element_blank()) +
  facet_wrap(~metabolite, ncol=5) +
  theme(legend.position="bottom")


ggsave("KX_Results_4_Conc_Stab_All.pdf", w=18, h=18)







#####################################################################################################################
################################################ Statistical Testing ################################################
#####################################################################################################################
 
 
 
rm(c_hi_long_Base,c_hi_long_BOTH,c_hi_long_FSBPase,c_hi_long_TKT)
rm(c_lo_long_Base,c_lo_long_BOTH,c_lo_long_TKT,c_lo_long_FSBPase)
rm(conc_hi_Base,conc_hi_BOTH,conc_hi_FSBPase,conc_hi_TKT)
rm(conc_lo_Base,conc_lo_BOTH,conc_lo_FSBPase,conc_lo_TKT)
rm(perc_lo_Base,perc_lo_BOTH,perc_lo_FSBPase,perc_lo_TKT)
rm(perc_hi_Base,perc_hi_BOTH,perc_hi_FSBPase,perc_hi_TKT)
rm(c_long_Base,c_long_BOTH,c_long_TKT, c_long_FSBPase)
rm(perc_data_Base,perc_data_BOTH,perc_data_TKT,perc_data_FSBPase)


Metabolite_List = unique(c_long_All$metabolite)
Nr_metabolites = length(Metabolite_List)

Nr_Models = length(unique(c_long_All$variant))
Nr_States = length(unique(c_long_All$stability))

Metabolite = Metabolite_List
KS_Test_Data_D = data.frame(Metabolite)
KS_Test_Data_P = data.frame(Metabolite)


KS_Base_D = c()
KS_FSBPase_D = c()
KS_TKT_D = c()
KS_Both_D = c()
KS_TKTBoth_D = c()
KS_FSBPaseBoth_D = c()
KS_FSBPaseTKT_D = c()

KS_Base_P = c()
KS_FSBPase_P = c()
KS_TKT_P = c()
KS_Both_P = c()
KS_TKTBoth_P = c()
KS_FSBPaseBoth_P = c()
KS_FSBPaseTKT_P = c()

for (i in 1:Nr_metabolites){
  # always the base model
  x = filter(c_long_All, metabolite == Metabolite_List[i], stability == "stable", variant == "4_Base")$concentration

  # FSBPase model
  y_FSBPase = filter(c_long_All, metabolite == Metabolite_List[i], stability == "stable", variant == "3_FSBPase")$concentration
  KS_Test_Result_FSBPase = ks.test(x, y_FSBPase, alternative="two.sided")
  KS_FSBPase_D[i] = KS_Test_Result_FSBPase$statistic
  KS_FSBPase_P[i] = KS_Test_Result_FSBPase$p.value

  # TKT
  y_TKT = filter(c_long_All, metabolite == Metabolite_List[i], stability == "stable", variant == "2_TKT")$concentration
  KS_Test_Result_TKT = ks.test(x, y_TKT, alternative="two.sided")
  KS_TKT_D[i] = KS_Test_Result_TKT$statistic
  KS_TKT_P[i] = KS_Test_Result_TKT$p.value

  # Both
  y_Both = filter(c_long_All, metabolite == Metabolite_List[i], stability == "stable", variant == "1_BOTH")$concentration
  KS_Test_Result_Both = ks.test(x, y_Both, alternative="two.sided")
  KS_Both_D[i] = KS_Test_Result_Both$statistic
  KS_Both_P[i] = KS_Test_Result_Both$p.value

  # TKT vs Both
  KS_Test_Result_TKTBoth = ks.test(y_TKT, y_Both, alternative="two.sided")
  KS_TKTBoth_D[i] = KS_Test_Result_TKTBoth$statistic
  KS_TKTBoth_P[i] = KS_Test_Result_TKTBoth$p.value

  # FSBPase vs Both
  KS_Test_Result_FSBPaseBoth = ks.test(y_FSBPase, y_Both, alternative="two.sided")
  KS_FSBPaseBoth_D[i] = KS_Test_Result_FSBPaseBoth$statistic
  KS_FSBPaseBoth_P[i] = KS_Test_Result_FSBPaseBoth$p.value

  # FSBPase vs TKT
  KS_Test_Result_FSBPaseTKT = ks.test(y_FSBPase, y_TKT, alternative="two.sided")
  KS_FSBPaseTKT_D[i] = KS_Test_Result_FSBPaseTKT$statistic
  KS_FSBPaseTKT_P[i] = KS_Test_Result_FSBPaseTKT$p.value

}


KS_Test_Data_D$BaseVsFSBPase = KS_FSBPase_D
KS_Test_Data_D$BaseVsTKT = KS_TKT_D
KS_Test_Data_D$BaseVsBoth = KS_Both_D
KS_Test_Data_D$TKTVsBoth = KS_TKTBoth_D
KS_Test_Data_D$FSBPasesVsBoth = KS_FSBPaseBoth_D
KS_Test_Data_D$FSBPaseVsTKT = KS_FSBPaseTKT_D


KS_Test_Data_P$BaseVsFBPase = KS_FSBPase_P
KS_Test_Data_P$BaseVsTKT = KS_TKT_P
KS_Test_Data_P$BaseVsBoth = KS_Both_P
KS_Test_Data_P$TKTVsBoth = KS_TKTBoth_P
KS_Test_Data_P$FSBPasesVsBoth = KS_FSBPaseBoth_P
KS_Test_Data_P$FSBPaseVsTKT = KS_FSBPaseTKT_P






KS_Test_Data_P_long = reshape2::melt(KS_Test_Data_P)
colnames(KS_Test_Data_P_long) = c("Metabolite", "Comparison","Pvalue")

Mets_To_Remove = c("O2","CO2_cax","CO2_cyt","COAPool","FUMPool","PPool")
for(i in 1:length(Mets_To_Remove)){
  KS_Test_Data_P_long = subset(KS_Test_Data_P_long, Metabolite!=Mets_To_Remove[i])
}


KS_Test_Data_P_long$Significant = rep("TRUE", length(KS_Test_Data_P_long$Metabolite))


for (i in 1:length(KS_Test_Data_P_long$Metabolite)){
  if(KS_Test_Data_P_long$Pvalue[i] > 0.05) {
    KS_Test_Data_P_long$Significant[i] = "FALSE"
  }
}

KS_Test_Data_P_long$Pvalue = round(KS_Test_Data_P_long$Pvalue, 3)



gp = ggplot(KS_Test_Data_P_long, aes(x=Comparison, y=Metabolite, fill=Significant))
gp = gp + geom_tile(colour="white", size=0.5)
gp = gp + geom_text(aes(label = Pvalue), size = 3)
gp = gp + theme_bw()
#gp = gp + coord_fixed()
gp = gp + scale_x_discrete(position="top")
gp
ggsave("KX_Results_5_KS-Test_Stable.pdf", w=10, h=10)
