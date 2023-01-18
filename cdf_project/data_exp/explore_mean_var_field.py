#%% Clean Memory and import pakages
for i in list(globals().keys()):
    if(i[0] != '_'):
        exec('del {}'.format(i))

from cdf_project.config import data_dir
from cdf_project.setup_testSet import setup_dict

import pickle
import numpy as np

#%%  Dataset paths
raw_data_dir = data_dir / 'raw'
field_data_dir = raw_data = data_dir / 'intermediate'
sol_data_dir = data_dir / 'processed'
point_data_dir = data_dir / 'models'






#mc_runs_array=[30, 20, 950]#setup_dict['monte_carlo_runs']
#mc_start_array =[0, 30, 50]
mc_runs_array=[500,500]#setup_dict['monte_carlo_runs']
mc_start_array =[0,500]


input_mesh = []
#%%  set fields file name
seeds_array_infiled = []
filed_hash_name= setup_dict['filed_hash_name']
filed_dir = f'fileds_{filed_hash_name}'



#%% merge fields

nxy = setup_dict['number_of_cells']
dims = setup_dict['dimensions']

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



    mean_field+=np.sum(fields, axis=0)
    variance_field+=np.sum(fields**2, axis=0)
    del fields


#%% merge fields
del fields_dict

N = len(seeds_array_infiled)
variance_sols = np.sqrt((variance_field - mean_field**2 / N)/(N-1))
mean_sols = mean_field / N



# %%

mean_matrix = mean_sols.reshape((nxy, nxy))
variance_matrix = variance_sols.reshape((nxy, nxy))

# %%
from explore_utils_fileds import show_mean_varaince_2d_3d
from cdf_project.data_etl.etl_utils import array_from_mesh
xy = array_from_mesh(dims, input_mesh)

fig = show_mean_varaince_2d_3d(mean_matrix, variance_matrix, xy, BW= False)



specifc=setup_dict['eq_formula2']
title_txt = f'All fields, Over:{N} \n {specifc}'  #
fig.suptitle(title_txt)
fig.tight_layout()
fig.show()


# %%


from cdf_project.config import reportings_dir

folder_name = setup_dict['id']
report_folder = reportings_dir/folder_name
try:
    report_folder.mkdir(parents=True, exist_ok=False)
except FileExistsError:
    print(f"Folder {report_folder} is already there")
else:
    print(f"Folder {report_folder} was created")


file_name = f'all_fileds_'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh']
fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")

# %%

# %%
