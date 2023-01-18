# -*- coding: utf-8 -*-
"""
Version 1, Exploration UTILS
Plots for showing fields or single solutions images and movies.
Plot by saving on file.
@author: enzo
"""
#%% import pakages
import matplotlib.pyplot as plt

#%% Set style
import matplotlib as mpl
axtickfsize = 16
labelfsize = 20
legfsize = labelfsize - 2
txtfsize = labelfsize - 2
lwidth = 3
markersize = 10
markeredgewidth = 0.1
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['axes.labelsize'] = labelfsize
mpl.rcParams['xtick.labelsize'] = axtickfsize
mpl.rcParams['ytick.labelsize'] = axtickfsize
mpl.rcParams['font.size'] = txtfsize
mpl.rcParams["figure.titlesize"] = 26
mpl.rcParams["figure.titleweight"] = 'regular'
mpl.rcParams['legend.fontsize'] = legfsize
mpl.rcParams['lines.linewidth'] = lwidth
mpl.rcParams['lines.markersize'] = markersize
mpl.rcParams['lines.markeredgewidth'] = markeredgewidth



#%% define

def min_max(x_seq, maxs, mins, title):
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))#
    
    ax.plot(x_seq, maxs, label='max')
    ax.plot(x_seq, mins, '.', label='min')
    ax.legend(loc="best")
    ax.grid(linestyle = '--', linewidth = 0.3)
    ax.set_xticks([int(j) for j in range(0,len(mins),int(len(mins)/20))])
    ax.set_title(title)
    #fig.show()
    return fig


#%% 7. def CDF generator

def C_p(sample):
    import numpy as np
    # calculate parameters
    sample_mean = np.mean(sample)
    sample_std = np.std(sample)
    print('Mean=%.3f, Standard Deviation=%.3f' % (sample_mean, sample_std))
    # calculate 
    data_sorted = np.sort(sample)
    # calculate the proportional values of samples
    p = 1. * np.arange(len(sample)) / (len(sample) - 1)
    return data_sorted, p

def cdf_1(center_values, x_label, y_label):
    # plot the sorted data:
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))#plt.figure(figsize=(10, 5))

    data_sorted, p = C_p(center_values) #n
    ax.plot(data_sorted, p, label=f'{x_label}; {y_label}')

    ax.set_xlabel('$h$') #SCALA 0- 1 E NO 0.2 - 0.24, PRENDERE IL C MAX OVVERO A 0.5
    ax.set_ylabel('$p$')
    #ax.set_title(f"{setup_dict['eq_formula']} \n {setup_dict['eq_formula2']}")

    Cmin = center_values.min()
    Cmax = center_values.max()
    print(f'h_min={Cmin}, h_max={Cmax}')
    left, right = ax.get_xlim()  # return the current xlim
    ax.set_xlim((min(left, Cmin), max(right, Cmax)))   # set the xlim to left, right
    
    import numpy as np
    ax.set_xticks(np.arange(Cmin, Cmax, step=(Cmax-Cmin)/5))
    ax.legend(loc="best")
    return fig




def cdf_2(center_values, x_label, y_label):

    fig = plt.figure(figsize=(16, 9)) # 10,5

    # plot the sorted data:

    ax1 = fig.add_subplot(121)
    data_sorted, p = C_p(center_values) #n
    ax1.plot(data_sorted, p, label=f'{x_label}; {y_label}')
    ax1.set_xlabel('$h$') #SCALA 0- 1 E NO 0.2 - 0.24, PRENDERE IL C MAX OVVERO A 0.5
    ax1.set_ylabel('$p$')
    ax1.set_title(f'relative')
    ax1.set_ylim(0,1)
    ax1.legend()

    ax2 = fig.add_subplot(122)
    data_sorted, p = C_p(center_values) #n
    Cmin = data_sorted.min()
    Cmax = data_sorted.max()
    data_sorted_N=(data_sorted-Cmin)/(Cmax - Cmin)
    ax2.plot(data_sorted_N, p)

    ax2.set_xlabel('$h$') #SCALA 0- 1 E NO 0.2 - 0.24, PRENDERE IL C MAX OVVERO A 0.5
    ax2.set_ylabel('$p$')
    ax2.set_title(f'normalized')
    ax2.set_ylim(0,1)
    import numpy as np
    ax2.set_xticks(np.arange(0, 1.01, step=1/10))
    ax2.legend()
    return fig


def cumulative(probe, mc_runs):
    fig, axs = plt.subplots(2, 2, figsize=(16, 9)) #10,10
    import numpy as np

    #mean
    cum_summ_dr = np.cumsum(probe)
    N = np.arange(2,mc_runs+1)
    N1 = np.arange(1,mc_runs)
    mean_csmdr = cum_summ_dr / N
    axs[0, 0].plot(N, mean_csmdr, 'tab:gray')
    axs[0, 0].set_title('Mean')

    #variance - second order
    var=probe- mean_csmdr
    var_res =  np.cumsum(var**2)/ N1
    std_res = np.sqrt(var_res)
    axs[0, 1].plot(N, var_res, 'tab:red')
    axs[0, 1].set_title('Variance')

    #skewness - third order
    skew_res =  np.cumsum((var/std_res)**3)/ N
    axs[1, 0].plot(N, skew_res, 'tab:green')
    axs[1, 0].set_title('Skewness')


    #kurtosis - 4th-order
    kurt_res =  np.cumsum((var/std_res)**4)/ N - 3
    axs[1, 1].plot(N, kurt_res, 'tab:orange')
    axs[1, 1].set_title('Kurtosis')


    return fig


















