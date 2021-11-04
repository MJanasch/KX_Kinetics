## Select the FCCs for a certain reaction from each sampling output file and concatenate them into one file
# This allows to calculate the median and MAD for the FCC for each reaction over all sampling results
# Needs to be looped over the number of reactions, here: 53

import sys, os
from glob import glob
import pandas as pd

Rxn_Nr = int(sys.argv[1])+2
Out_File = sys.argv[2]

def read_col(fn):
	return pd.read_csv(fn, sep="\t", usecols=[0, 1, 2, Rxn_Nr], header=None)

files = glob('../FOLDER_CONTAINING_FCCs/FCCs_DATE_fMCS_*.tab')

#files = ['/hdd/markus/KX/KX_parameter_sampling/KX_parameter_sampling_results_211013_TKT/FCCs_211013_fMCS_1.tab', '/hdd/markus/KX/KX_parameter_sampling/KX_parameter_sampling_results_211013_TKT/FCCs_211013_fMCS_2.tab']

big_df = pd.concat([read_col(fn) for fn in files], axis=0)






big_df.to_csv(Out_File, sep ='\t')