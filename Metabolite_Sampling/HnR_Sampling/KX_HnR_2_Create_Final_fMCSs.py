## Random Sampling of the Hit-and-Run results of the 87° HnR sets
#° depending on the number of HnR out files, which depends on the number of 
# initial concentration sets from the metabolite variability analysis, about 2x Nr of metabolites

# Import libraries
import sys, os
import numpy as np
import time
import math
from scipy import stats
import random

# Start timer
start_time = time.time()


## Read in datasets and append into one lHe file (should be around 400K)
# Extract header from first HnR sampling file
Conc_File_Header_name = "HnR_Sampling_File_1.tsv"
Conc_File_Header_content= open(Conc_File_Header_name, "r")
header = Conc_File_Header_content.readline().rstrip()



## Loop through HnR concentration sets, discard header, save concentrations
Nr_Starting_Set = range(1,87)

Total_Content = ""
Nr_Datasets = 0
for Nr_Set in Nr_Starting_Set:
	print(Nr_Set)
	Conc_File_name = "HnR_Sampling_File_"+str(Nr_Set)+".tsv"
	Conc_File_content= open(Conc_File_name, "r")
	Content = ""
	for line in Conc_File_content:
		line=line.rstrip()
		if line[0] != "R":
			Content = Content + line + "\n"
			Nr_Datasets += 1
	Total_Content = Total_Content + Content
	#print(Total_Content)
	print(Nr_Datasets)
	print("--- %s seconds ---" % (time.time() - start_time))

## Create 10000 random numbers between 0 and the total Number of Datasets
# Repeat 5 times to check for bias
rn_nrs = []
count = 0
while count < 10000:
	rn_nr = random.randint(0,Nr_Datasets)
	if not rn_nr in rn_nrs:
		rn_nrs.append(rn_nr)
		count += 1
#print(rn_nrs)


rn_nrs_2 = []
count_2 = 0
while count_2 < 10000:
	rn_nr_2 = random.randint(0,Nr_Datasets)
	if not rn_nr_2 in rn_nrs_2:
		rn_nrs_2.append(rn_nr_2)
		count_2 += 1
#print(rn_nrs_2)

rn_nrs_3 = []
count_3 = 0
while count_3 < 10000:
	rn_nr_3 = random.randint(0,Nr_Datasets)
	if not rn_nr_3 in rn_nrs_3:
		rn_nrs_3.append(rn_nr_3)
		count_3 += 1
#print(rn_nrs_3)


rn_nrs_4 = []
count_4 = 0
while count_4 < 10000:
	rn_nr_4 = random.randint(0,Nr_Datasets)
	if not rn_nr_4 in rn_nrs_4:
		rn_nrs_4.append(rn_nr_4)
		count_4 += 1
#print(rn_nrs_4)

rn_nrs_5 = []
count_5 = 0
while count_5 < 10000:
	rn_nr_5 = random.randint(0,Nr_Datasets)
	if not rn_nr_5 in rn_nrs_5:
		rn_nrs_5.append(rn_nr_5)
		count_5 += 1
#print(rn_nrs_5)

## Take the fMCSs corresponding to the random numbers in rn_nrs
fMCSs = list(Total_Content.split("\n"))
Final_fMCSs = ""
for nr in rn_nrs:
	Final_fMCSs = Final_fMCSs + fMCSs[nr] + "\n"

Final_fMCSs_2 = ""
for nr_2 in rn_nrs_2:
	Final_fMCSs_2 = Final_fMCSs_2 + fMCSs[nr_2] + "\n"

Final_fMCSs_3 = ""
for nr_3 in rn_nrs_3:
	Final_fMCSs_3 = Final_fMCSs_3 + fMCSs[nr_3] + "\n"

Final_fMCSs_4 = ""
for nr_4 in rn_nrs_4:
	Final_fMCSs_4 = Final_fMCSs_4 + fMCSs[nr_4] + "\n"

Final_fMCSs_5 = ""
for nr_5 in rn_nrs_5:
	Final_fMCSs_5 = Final_fMCSs_5 + fMCSs[nr_5] + "\n"



## Output
outfile_name = "HnR_SubSampled_Output.txt" # file for further analysis in matlab in .tsv file format
outfile = open(outfile_name, "w")
outfile.write(Final_fMCSs)
outfile.close()

outfile_name_2 = "HnR_SubSampled_Output_2.txt" # file for further analysis in matlab in .tsv file format
outfile_2 = open(outfile_name_2, "w")
outfile_2.write(Final_fMCSs_2)
outfile_2.close()

outfile_name_3 = "HnR_SubSampled_Output_3.txt" # file for further analysis in matlab in .tsv file format
outfile_3 = open(outfile_name_3, "w")
outfile_3.write(Final_fMCSs_3)
outfile_3.close()

outfile_name_4 = "HnR_SubSampled_Output_4.txt" # file for further analysis in matlab in .tsv file format
outfile_4 = open(outfile_name_4, "w")
outfile_4.write(Final_fMCSs_4)
outfile_4.close()

outfile_name_5 = "HnR_SubSampled_Output_5.txt" # file for further analysis in matlab in .tsv file format
outfile_5 = open(outfile_name_5, "w")
outfile_5.write(Final_fMCSs_5)
outfile_5.close()
