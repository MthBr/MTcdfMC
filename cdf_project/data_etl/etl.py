# -*- coding: utf-8 -*-
"""
Version 1, Extraction Transformation Loading
Read setup, generete and write / or read fileds
Genereate and write single (1) solution
@author: enzo
"""
#%% Clean Memory and import pakages
for i in list(globals().keys()):
    if(i[0] != '_'):
        exec('del {}'.format(i))

from cdf_project.config import data_dir
from cdf_project.setup_testSet import setup_dict
from cdf_project.data_etl.etl_utils import gen_fields

import pickle

#%%  Dataset paths
raw_data_dir = data_dir / 'raw'
field_data_dir = raw_data = data_dir / 'intermediate'
single_run_data_dir = data_dir / 'intermediate'


#%%  set fields file name

dimens = setup_dict['dimensions']
crrl_value=setup_dict['correlation_coeff'] 
n = setup_dict['number_of_cells']
mc_runs=500#setup_dict['monte_carlo_runs']
mc_start = 500
i_leng = setup_dict['length']
srt = setup_dict['start']
is_log_normal= setup_dict['log_normal']

filed_hash_name= setup_dict['filed_hash_name']


filed_dir = f'fileds_{filed_hash_name}'
(field_data_dir / filed_dir).mkdir(exist_ok=True) #parents=True, 


field_file_name=f'fileds_normal_mc{mc_runs}_s{mc_start}_{filed_hash_name}.plk'

fields_file = field_data_dir / filed_dir / field_file_name


#%%  if field exists read evals genererate
if fields_file.is_file():
    print('Exist, loading')
    with open(fields_file, 'rb') as f:
        fields = pickle.load(f)
else:
    print('Not exist, generating')
    fields = gen_fields(dims= dimens, number_of_cells=n, input_dict=setup_dict, ens_num=mc_runs, enstart= mc_start, log_normal=is_log_normal, interval_len = i_leng, interval_start= srt)
    with open(fields_file, 'wb') as f:
        pickle.dump(fields, f, protocol=4)#TODO higer protcol, ...





#%%  Set solution number, and solution name

#with open(fields_file, 'wb') as f:
#    pickle.dump(fields, f, protocol=4)#TODO higer protcol, ...


#%%  Generate and write solution (partial savings, with intervals ?)







