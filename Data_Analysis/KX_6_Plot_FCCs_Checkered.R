library(tidyverse)
library(data.table)
# Load data
#example_data = read_tsv("cbb_example_data.tab", na="NA")

Base_data = as.data.frame(fread("KX_Results_3_FCCs_Med_MAD_Base.csv"))
Base_data = subset(Base_data, select = -(V1))
Base_data$Model = rep("Base", nrow(Base_data))
Base_data$Has_TKT_Reg = rep("No", nrow(Base_data))
Base_data$Has_FSBPase_Reg = rep("No", nrow(Base_data))

FSBPase_data = as.data.frame(fread("KX_Results_3_FCCs_Med_MAD_FSBPase.csv"))
FSBPase_data = subset(FSBPase_data, select = -(V1))
FSBPase_data$Model = rep("FSBPase", nrow(FSBPase_data))
FSBPase_data$Has_TKT_Reg = rep("No", nrow(FSBPase_data))
FSBPase_data$Has_FSBPase_Reg = rep("Yes", nrow(FSBPase_data))

TKT_data = as.data.frame(fread("KX_Results_3_FCCs_Med_MAD_TKT.csv"))
TKT_data = subset(TKT_data, select = -(V1))
TKT_data$Model = rep("TKT", nrow(TKT_data))
TKT_data$Has_TKT_Reg = rep("Yes", nrow(TKT_data))
TKT_data$Has_FSBPase_Reg = rep("No", nrow(TKT_data))

Both_data = as.data.frame(fread("KX_Results_3_FCCs_Med_MAD_BOTH.csv"))
Both_data = subset(Both_data, select = -(V1))
Both_data$Model = rep("Both", nrow(Both_data))
Both_data$Has_TKT_Reg = rep("Yes", nrow(Both_data))
Both_data$Has_FSBPase_Reg = rep("Yes", nrow(Both_data))

Data_long_All = rbind(Base_data,FSBPase_data,TKT_data,Both_data)



## Sort reactions and group them according to the pathway
Data_long_All[Data_long_All == "RuBisC"] = "01_RuBisC"
Data_long_All[Data_long_All == "PGK"] = "02_PGK"
Data_long_All[Data_long_All == "GAPD"] = "03_GAPD"
Data_long_All[Data_long_All == "TPI"] = "04_TPI"
Data_long_All[Data_long_All == "ALD"] = "05_ALD"
Data_long_All[Data_long_All == "FBPase"] = "06_FBPase"
Data_long_All[Data_long_All == "TKT1"] = "07_TKT1"
Data_long_All[Data_long_All == "TKT2"] = "08_TKT2"
Data_long_All[Data_long_All == "TAL"] = "09_TAL"
Data_long_All[Data_long_All == "FBA"] = "10_FBA"
Data_long_All[Data_long_All == "SBPase"] = "11_SBPase"
Data_long_All[Data_long_All == "RPE"] = "12_RPE"
Data_long_All[Data_long_All == "RPI"] = "13_RPI"
Data_long_All[Data_long_All == "PRK"] = "14_PRK"
Data_long_All[Data_long_All == "PGM"] = "15_PGM"
Data_long_All[Data_long_All == "ENO"] = "16_ENO"
Data_long_All[Data_long_All == "PYK"] = "17_PYK"
Data_long_All[Data_long_All == "PDH"] = "18_PDH"
Data_long_All[Data_long_All == "ME"] = "19_ME"
Data_long_All[Data_long_All == "PPK"] = "20_PPK"

Data_long_All[Data_long_All == "CS"] = "21_CS"
Data_long_All[Data_long_All == "FUMC"] = "22_FUMC"
Data_long_All[Data_long_All == "ICD"] = "23_ICD"
Data_long_All[Data_long_All == "MDH"] = "24_MDH"
Data_long_All[Data_long_All == "ACO"] = "25_ACO"

Data_long_All[Data_long_All == "ATPSyn"] = "26_ATPSyn"
Data_long_All[Data_long_All == "NADPase"] = "27_NADPase"
Data_long_All[Data_long_All == "Supply_P_i"] = "28_Supply_Pi"
Data_long_All[Data_long_All == "RuBisO"] = "29_RuBisO"
Data_long_All[Data_long_All == "PGP"] = "30_PGP"
Data_long_All[Data_long_All == "GLCO"] = "31_GLCO"
Data_long_All[Data_long_All == "PGI"] = "32_PGI"
Data_long_All[Data_long_All == "ZWF"] = "33_ZWF"
Data_long_All[Data_long_All == "DEVB"] = "34_DEVB"
Data_long_All[Data_long_All == "GND"] = "35_GND"

Data_long_All[Data_long_All == "PRSA"] = "36_PRSA"
Data_long_All[Data_long_All == "ADK"] = "37_ADK"

Data_long_All[Data_long_All == "XFPK1"] = "38_XFPK1"
Data_long_All[Data_long_All == "XFPK2"] = "39_XFPK2"
Data_long_All[Data_long_All == "PTA"] = "40_PTA"


Rxn_To_Remove = c("RuBisO","PGP","GLCO","Sink_E4P","Sink_ACCOA","Sink_G6P","Sink_GLX","Sink_H2O2","Sink_O2G","Sink_OAA","Sink_P3G","Sink_PRPP","Sink_PYR","Supply_COA","Supply_FUM","GND","PGI","ZWF","DEVB","PNT","PPK","XFPK1","XFPK2","PTA","ACO","PRSA","ADK")
Data_long_All_Mod = Data_long_All
for(i in 1:length(Rxn_To_Remove)){
  Data_long_All_Mod = subset(Data_long_All_Mod, Effector!=Rxn_To_Remove[i])
}

for(i in 1:length(Rxn_To_Remove)){
  Data_long_All_Mod = subset(Data_long_All_Mod, Target!=Rxn_To_Remove[i])
}


# Define organism color scale
models = c("Base", "FSBPase", "TKT", "Both")
organcols = c("#762a83", "#9970ab","#5aae61","#1b7837")
names(organcols) = models






gp_supp = ggplot(Data_long_All, aes(x=Has_FSBPase_Reg, y=reorder(Has_TKT_Reg,desc(Has_TKT_Reg)), fill=Median_FCC, alpha=MAD))
gp_supp = gp_supp + geom_tile()
gp_supp = gp_supp + facet_grid(Effector~Target,switch="both")
gp_supp = gp_supp + scale_fill_gradientn(colours=c("#fdae61","#ffffbf","#2c7bb6"), values=rescale(c(min(Data_long_All$Median_FCC),0,max(Data_long_All$Median_FCC)), c(0,1)))
gp_supp = gp_supp + scale_alpha_continuous(range=c(1,0.5))
gp_supp = gp_supp + theme_bw()
gp_supp = gp_supp + theme(
  strip.background = element_blank(),
  aspect.ratio = 1,
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  strip.text.y.left = element_text(angle=0, hjust=1, vjust=0.5, size = 12),
  strip.text.x.bottom= element_text(angle=90, hjust=1, vjust=1, size = 12),
  panel.spacing = unit(0.1, "lines"),
  panel.border = element_rect(color="#d9d9d9")
)
#gp_supp
# Save the plot
ggsave("KX_Results_6_All_FCCs_Supplement.pdf", w=18, h=18)



# Make small plot for sekected Reactions
gp = ggplot(Data_long_All_Mod, aes(x=Has_FSBPase_Reg, y=reorder(Has_TKT_Reg,desc(Has_TKT_Reg)), fill=Median_FCC, alpha=MAD))
gp = gp + geom_tile()
gp = gp + facet_grid(Effector~Target,switch="both")
gp = gp + scale_fill_gradientn(colours=c("#fdae61","#ffffbf","#2c7bb6"), values=rescale(c(min(Data_long_All$Median_FCC),0,max(Data_long_All$Median_FCC)), c(0,1)))
gp = gp + scale_alpha_continuous(range=c(1,0.5))
gp = gp + theme_bw()
gp = gp + theme(
  strip.background = element_blank(),
  aspect.ratio = 1,
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  strip.text.y.left = element_text(angle=0, hjust=1, vjust=0.5, size = 12),
  strip.text.x.bottom= element_text(angle=90, hjust=1, vjust=1, size = 12),
  panel.spacing = unit(0.1, "lines"),
  panel.border = element_rect(color="#d9d9d9")
)
#ggsave("quad_plot_example.pdf", w=9, h=4.5)
