
#%% Import
import numpy as np


#%% 7. def CDF generator

def C_p(allMC_phi_tt, m, n):
    sample = allMC_phi_tt[:, m, n]  # tempo x spazio
    # calculate parameters
    sample_mean = np.mean(sample)
    sample_std = np.std(sample)
    print('time = %.d   Mean=%.3f, Standard Deviation=%.3f' % (m, sample_mean, sample_std))
    # calculate 
    data_sorted = np.sort(sample)

    # calculate the proportional values of samples
    p = 1. * np.arange(len(sample)) / (len(sample) - 1)
    return data_sorted, p


