# -*- coding: utf-8 -*-
"""
Version 1, Exploration
Read setup, plot fileds
Read and plot single (1) solution
@author: enzo
"""
#%% Clean Memory and import pakages
for i in list(globals().keys()):
    if(i[0] != '_'):
        exec('del {}'.format(i))

#%%
from cdf_project.config import data_dir
from cdf_project.setup_testSet import setup_dict
import matplotlib.pyplot as plt
import pickle

#%%  Dataset paths
field_data_dir = raw_data = data_dir / 'intermediate'
#single_run_data_dir = data_dir / 'intermediate'


#%%  set fields file name and soltion file name
dims = setup_dict['dimensions']
crrl_value=setup_dict['correlation_coeff'] 
crrl_text=setup_dict['correl_text'] 
n = setup_dict['number_of_cells']

len = setup_dict['length']
srt = setup_dict['start']


mc_runs=500#setup_dict['monte_carlo_runs']
mc_start = 500


filed_hash_name= setup_dict['filed_hash_name']
filed_dir = f'fileds_{filed_hash_name}'
field_file_name=f'fileds_normal_mc{mc_runs}_s{mc_start}_{filed_hash_name}.plk'
fields_file = field_data_dir / filed_dir / field_file_name


#%%  Read fileds and solution
if fields_file.is_file():
    print('Exist, loading')
    with open(fields_file, 'rb') as f:
        fields_dict = pickle.load(f)
else:
    print('Not exists!!!!')


#%% 
fields_array = fields_dict['fields_array']
input_mesh = fields_dict['mesh']
seeds_array = fields_dict['seeds_array']
del fields_dict
#%% 
from cdf_project.data_etl.etl_utils import gen_uniform_delta
delta_space = gen_uniform_delta(n)  
setup_dict['grid_spacing'] = delta_space


#%% Plot selected fields
indx= 0

field= fields_array[indx]


#from cdf_project.data_exp.explore_utils_fileds import plot4fileds
#fig, ax = plot4fileds(field, dims)


#%% Calculate selected field slution
from cdf_project.feature_eng.solver4stationary_util import single_solve

import json

name_cnf = 'AMG_CLASSICAL_PMIS'
file_conf = name_cnf+'.json'
f = open(f"{data_dir/'configs'/file_conf}")
cfg_dict = json.load(f)
setup_dict['delta_type'] = '2-Cubic' #
final_phi  = single_solve(setup_dict, input_mesh, field.flatten(), cfg_dict)  #False  True

#%% 
import numpy as np
maxs = np.max(final_phi)
mins = np.min(final_phi)
print("Values bigger than 0 =", maxs[maxs>=0])
print("Their indices are ", np.nonzero(maxs >= 0))

nans = np.count_nonzero(np.isnan(final_phi))
print("Values of nan", nans)




#%%  Plot selected solution
from cdf_project.data_etl.etl_utils import array_from_mesh
xy=array_from_mesh(2, input_mesh)


values_matrix = final_phi.reshape((n, n))

from cdf_project.data_exp.explore_utils_sols import show_sinlge_2d_sol_filed, show_sinlge_2d_sol_filed_in3d

formula=setup_dict['eq_formula']
specifc=setup_dict['eq_formula2']
title_txt = f'{formula} \n {specifc}'  #
fig = show_sinlge_2d_sol_filed_in3d(field, values_matrix, xy, BW=False)
fig.suptitle(title_txt)
fig.tight_layout()




#%%  Plot selected solution

from cdf_project.config import reportings_dir
folder_name = setup_dict['id']
report_folder = reportings_dir/folder_name
try:
    report_folder.mkdir(parents=True, exist_ok=False)
except FileExistsError:
    print(f"Folder {report_folder} is already there")
else:
    print(f"Folder {report_folder} was created")

# %%
file_name = f'{indx+mc_start}_s{seeds_array[indx]}_{indx}'+f'{name_cnf}'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq'] + '_3d'

fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")



# %%

# %%
