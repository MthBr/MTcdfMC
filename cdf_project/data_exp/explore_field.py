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


from cdf_project.config import data_dir
from cdf_project.setup_testSet import setup_dict
import matplotlib.pyplot as plt
import pickle

#%%  Dataset paths
field_data_dir = raw_data = data_dir / 'intermediate'
single_run_data_dir = data_dir / 'intermediate'


#%%  set fields file name and soltion file name

#%%  set fields file name

dimens = setup_dict['dimensions']
crrl_value=setup_dict['correlation_coeff'] 
crrl_text=setup_dict['correl_text'] 
n = setup_dict['number_of_cells']

mc_runs_array=[500]#setup_dict['monte_carlo_runs']
mc_start_array =[0]

seeds_array_infiled = []
filed_hash_name= setup_dict['filed_hash_name']
filed_dir = f'fileds_{filed_hash_name}'

#%% merge fields

nxy = setup_dict['number_of_cells']
dims = setup_dict['dimensions']

import numpy as np

mean_field= np.zeros((nxy,nxy))
variance_field=np.zeros((nxy,nxy))


#%% merge fields
for mc_runs, mc_start in zip(mc_runs_array,mc_start_array):
    #%%  set fields file name
    field_file_name=f'fileds_normal_mc{mc_runs}_s{mc_start}_{filed_hash_name}.plk'
    fields_file = field_data_dir / filed_dir / field_file_name
    #%%  if field exists read evals genererate
    if fields_file.is_file():
        print('Exist, loading')
        with open(fields_file, 'rb') as f:
            fields_dict = pickle.load(f)
    else:
        print('Not exists!!!!')

    #%%  if field exists read evals genererate
    input_mesh = fields_dict['mesh']
    seeds_array_infiled = np.append(seeds_array_infiled, fields_dict['seeds_array'])   
    #np.alltrue(input_mesh.cellCenters.value == input_mesh.cellCenters.value)
    fields = fields_dict['fields_array']

#fields, delta_space, d_vector_mesh
#%% Plot selected fields

from cdf_project.data_exp.explore_utils_fileds import plot25fileds
filed_list = fields[0 : 25]

fig, ax = plot25fileds(filed_list)



#%% Plot selected fields
filed_list = fields[0 : 4]


from cdf_project.data_exp.explore_utils_fileds import plot4fileds
fig, ax = plot4fileds(filed_list, dimens)














# %%
