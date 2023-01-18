#%% Clean Memory and import pakages
for i in list(globals().keys()):
    if(i[0] != '_'):
        exec('del {}'.format(i))

from cdf_project.config import data_dir
from cdf_project.setup_testSet import setup_dict
from cdf_project.feature_eng.solver4stationary_util import ensamble_solutions

import pickle

#%%  Dataset paths
field_data_dir = raw_data = data_dir / 'intermediate'
sol_data_dir = data_dir / 'processed'


#%%  set fields file name

mc_runs=500#setup_dict['monte_carlo_runs']
mc_start = 500
filed_hash_name= setup_dict['filed_hash_name']


filed_dir = f'fileds_{filed_hash_name}'
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


fields_array = fields_dict['fields_array']
input_mesh = fields_dict['mesh']
seeds_array = fields_dict['seeds_array']
del fields_dict

from cdf_project.data_etl.etl_utils import gen_uniform_delta
n = setup_dict['number_of_cells']
delta_space = gen_uniform_delta(n)  
setup_dict['grid_spacing'] = delta_space


#%% Save or create
sol_hash_name = setup_dict['sol_hash_name']
sol_dir = f'sols_{sol_hash_name}'
(sol_data_dir / sol_dir).mkdir(exist_ok=True) #parents=True, 


solution_file_name=f'solution_mc{mc_runs}_s{mc_start}_{sol_hash_name}.plk'

solutions_file = sol_data_dir / sol_dir / solution_file_name


#%%  if solutions exists read evals genererate
if solutions_file.is_file():
    print('Exist, loading')
    with open(solutions_file, 'rb') as f:
        final_phi_arrays_dict = pickle.load(f)
else:
    print('Not exist, generating')
    final_phi_arrays_dict  = ensamble_solutions(setup_dict, input_mesh, fields_array)  #False  True
    final_phi_arrays_dict['seeds_array']=seeds_array
    with open(solutions_file, 'wb') as f:
        pickle.dump(final_phi_arrays_dict, f)




# %%
sols = final_phi_arrays_dict['final_solutions_array'] 
seeds_array = final_phi_arrays_dict['seeds_array']
del final_phi_arrays_dict


import random
lst_img_indx = random.sample(range(0, mc_runs), 5)
#lst_img_indx = range(0,mc_runs)
#lst_img_indx = range(0,15)
#indx=5


from cdf_project.data_etl.etl_utils import array_from_mesh
xy=array_from_mesh(2, input_mesh)


#%% SAVE

formula=setup_dict['eq_formula']
specifc=setup_dict['eq_formula2']
title_txt = f'{formula} \n {specifc}'  #

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

for indx in lst_img_indx:
    from cdf_project.data_exp.explore_utils_sols import show_sinlge_2d_sol_filed, show_sinlge_2d_sol_filed_in3d
    dimens = setup_dict['dimensions']
    nxy = setup_dict['number_of_cells']

    values_matrix = sols[indx].reshape((nxy, nxy))
    field=fields_array[indx]
    import numpy as np
    ind = [nxy//2,nxy//2]

    fig = show_sinlge_2d_sol_filed_in3d(field, values_matrix, xy, BW=True)
    fig.suptitle(title_txt)
    fig.tight_layout()
    file_name = f'{indx+mc_start}_s{seeds_array[indx]}_{indx}'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq'] + '_3dbw'
    
    import pylab as plt
    # for pubblication: Nature: Postscript, Vector EPS or PDF format
    #maximum 300 p.p.i.; and min 300dpi
    fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")
    # fig.savefig(report_folder/(file_name+'.pdf'), bbox_inches="tight")
    # fig.savefig(report_folder/(file_name+'.png'), bbox_inches="tight")
    # fig.savefig(report_folder/(file_name+'.eps'), format='eps', dpi=300)
    # fig.savefig(report_folder/(file_name+'.svg'), format='svg', dpi=300)
    plt.close(fig)

    fig = show_sinlge_2d_sol_filed_in3d(field, values_matrix, xy, BW=False)
    fig.suptitle(title_txt)
    fig.tight_layout()
    file_name = f'{indx+mc_start}_s{seeds_array[indx]}_{indx}'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq'] + '_3d'
    fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")
    plt.close(fig)




    fig, ax = show_sinlge_2d_sol_filed(field, values_matrix, xy, ind)
    fig.suptitle(title_txt)
    fig.tight_layout()
    file_name = f'{indx+mc_start}_s{seeds_array[indx]}_{indx}'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq']
    fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")
    plt.close(fig)



# %%
