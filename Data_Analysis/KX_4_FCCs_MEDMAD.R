#!/usr/bin/env Rscript



# 1 Load reaction header
# 2 Loop through reactions
  # 3 read in FCC files
  # 4 make FCCs close to 0 to 0
  # 5 calculate median and MAD
  # 6 put into dataframe
  # 7 concatenate
# give out csv file
# Whoop whoop, works!



# Read infile names from command line

header_file = "KX_Reaction_Header.txt" # Reaction header file


# Load data
library(data.table)
library(reshape2)

# The enzyme in the column influences the flux of reactions in rows (FCC)

# Load custom labels for reactions
custom_rxn_labels = scan(header_file, character(), quote = "")
length(custom_rxn_labels)

# Initialize dataframe
FCCS_File_Ind = paste("FCC_Concat_Reaction_1.tab",sep="",collapse="")
fccs_ind = as.data.frame(fread(FCCS_File_Ind))
fccs_ind[,4:ncol(fccs_ind)][abs(fccs_ind[,4:ncol(fccs_ind)]) < 1e-6] = 0 # Markus approves
fccs_ind = subset(fccs_ind, select = -(V1))
fccs_ind = fccs_ind[-c(1),]
colnames(fccs_ind) = c("Conc_set","Stable_set","Reaction",custom_rxn_labels[1])
fccs_ind$Reaction = custom_rxn_labels[fccs_ind$Reaction]
fccs_ind_med_mad_list = lapply(colnames(fccs_ind)[4:ncol(fccs_ind)], function(effector){
  # Calculate Median and MAD
  effector_median = aggregate(fccs_ind[,effector], list(Target = fccs_ind[,"Reaction"]), median)
  effector_mad = aggregate(fccs_ind[,effector], list(Target = fccs_ind[,"Reaction"]), mad)
  # The second column is the Median or the MAD
  colnames(effector_median)[2] = "Median_FCC"
  colnames(effector_mad)[2] = "MAD"
  # Add the Effector
  effector_median$Effector = effector
  effector_mad$Effector = effector
  #
  # Return the merged data frame
  merge(effector_median, effector_mad)
})

fccs_med_mad_ALL = as.data.frame(rbindlist(fccs_ind_med_mad_list))




for(i in 2:length(custom_rxn_labels)){
  FCCS_File_Ind = paste("FCC_Concat_Reaction_",i,".tab",sep="",collapse="")
  fccs_ind = as.data.frame(fread(FCCS_File_Ind))
  fccs_ind = subset(fccs_ind, select = -(V1))
  fccs_ind = fccs_ind[-c(1),]
  colnames(fccs_ind) = c("Conc_set","Stable_set","Reaction",custom_rxn_labels[i])
  fccs_ind[,4:ncol(fccs_ind)][abs(fccs_ind[,4:ncol(fccs_ind)]) < 1e-6] = 0 # Markus approves
  fccs_ind$Reaction = custom_rxn_labels[fccs_ind$Reaction]
  fccs_ind_med_mad_list = lapply(colnames(fccs_ind)[4:ncol(fccs_ind)], function(effector){
    # Calculate Median and MAD
    effector_median = aggregate(fccs_ind[,effector], list(Target = fccs_ind[,"Reaction"]), median)
    effector_mad = aggregate(fccs_ind[,effector], list(Target = fccs_ind[,"Reaction"]), mad)
    # The second column is the Median or the MAD
    colnames(effector_median)[2] = "Median_FCC"
    colnames(effector_mad)[2] = "MAD"
    # Add the Effector
    effector_median$Effector = effector
    effector_mad$Effector = effector
    #
    # Return the merged data frame
    merge(effector_median, effector_mad)
  })
  fccs_ind_med_mad = as.data.frame(rbindlist(fccs_ind_med_mad_list))
  
  fccs_med_mad_ALL = rbind(fccs_med_mad_ALL,fccs_ind_med_mad)
  
}

write.csv(fccs_med_mad_ALL,'KX_Results_3_FCCs_Med_MAD.csv')