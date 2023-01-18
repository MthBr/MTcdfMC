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


#%%  set MC variables

#mc_runs_array=[30, 20, 950]#setup_dict['monte_carlo_runs']
#mc_start_array =[0, 30, 50]
mc_runs_array=[500,500]#setup_dict['monte_carlo_runs']
mc_start_array =[0,500]


input_mesh = []
#%%  set fields file name
seeds_array_infiled = []
filed_hash_name= setup_dict['filed_hash_name']
filed_dir = f'fileds_{filed_hash_name}'

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

del fields_dict
#%% merge fields
#np.alltrue(seeds_array_infield0[-1] ==seeds_array_infield[0])
#np.append(seeds_array_infield0, seeds_array_infield)
#np.append(fields_array0, fields_array, axis=0)

nxy = setup_dict['number_of_cells']
dims = setup_dict['dimensions']

mean_sols= np.zeros(nxy**dims)
variance_sols=np.zeros(nxy**dims)



#%% Save or create
sol_hash_name = setup_dict['sol_hash_name']
sol_dir = f'sols_{sol_hash_name}'
maxs = []
mins =[]
center_values = []
seeds_array = []

for mc_runs, mc_start in zip(mc_runs_array,mc_start_array):
    solutions_file_name=f'solution_mc{mc_runs}_s{mc_start}_{sol_hash_name}.plk'
    solutions_file = sol_data_dir / sol_dir / solutions_file_name

    #%%  if solutions exists read evals genererate
    if solutions_file.is_file():
        print('Exist, loading')
        with open(solutions_file, 'rb') as f:
            solutions_dict = pickle.load(f)
    else:
        print('Not exist!!!!')

    seeds_array= np.append(seeds_array, solutions_dict['seeds_array'])   
    sols = solutions_dict['final_solutions_array'] 
    maxs =  np.append(maxs, np.max(sols, axis=1) )
    mins = np.append(mins, np.min(sols, axis=1))


    mean_sols+=np.sum(sols, axis=0)
    variance_sols+=np.sum(sols**2, axis=0)

    ##
    values_matrix = sols.reshape((len(sols),nxy, nxy))
    center_values = np.append(center_values, values_matrix[:, nxy//2, nxy//2])
    del solutions_dict
    del values_matrix

# mean_sols is sum
# variance_sols is sum**2
N = len(center_values)
assert (N == len(seeds_array_infiled))

variance_sols = np.sqrt((variance_sols - mean_sols**2 / N)/(N-1))
mean_sols = mean_sols / N


#%%  if sol exists read evals genererate
#sols = solutions_dict['final_solutions_array'] 
#seeds_array_insol= solutions_dict['seeds_array']

# Assert quality, merge, and extract central point!
# Plot max and min over solutions!!!!
#maxs = np.max(sols, axis=1)

print("Values bigger than 0 =", maxs[maxs>=0])
print("Their indices are ", np.nonzero(maxs >= 0))

#mins = np.min(sols, axis=1)
#mins






# %%

x_seq = np.arange(len(mins))

from util_stats_pdf import min_max
title=f"{setup_dict['eq_formula']} \n {setup_dict['eq_formula2']}"
fig = min_max(x_seq, maxs, mins, title)

fig.tight_layout()
from cdf_project.config import reportings_dir

file_name = f'minmax'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq']
folder_name = setup_dict['id']
report_folder = reportings_dir/folder_name
try:
    report_folder.mkdir(parents=True, exist_ok=False)
except FileExistsError:
    print(f"Folder {report_folder} is already there")
else:
    print(f"Folder {report_folder} was created")


# for pubblication: Nature: Postscript, Vector EPS or PDF format
#maximum 300 p.p.i.; and min 300dpi
fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")




# %%

mean_matrix = mean_sols.reshape((nxy, nxy))
variance_matrix = variance_sols.reshape((nxy, nxy))


from cdf_project.data_etl.etl_utils import array_from_mesh
xy = array_from_mesh(dims, input_mesh)


from matplotlib import ticker
#norm = LogNorm()

# from cdf_project.data_exp.explore_utils_fileds import show_mean_varaince_2d_3d
# fig = show_mean_varaince_2d_3d(mean_matrix, variance_matrix, xy, BW= False,\
#     loct={'mean':ticker.AutoLocator(), 'var':ticker.AutoLocator()})


#%%

from cdf_project.data_exp.explore_utils_fileds import show_mean_varaince_LOG_2d_3d
fig = show_mean_varaince_LOG_2d_3d(mean_matrix, variance_matrix, xy, BW= False)



#%%

formula=setup_dict['eq_formula']
specifc=setup_dict['eq_formula2']
title_txt = f'All Solutions Over:{N} \n {formula}  \n {specifc}'  #
fig.suptitle(title_txt)
fig.tight_layout()
fig.show()


file_name = f'all_solutionsMeanVar_'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh']
fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")


# %% Centered
def crop_center(matrix,perc=0.3):
    y,x = matrix.shape
    cropx= int(x*perc)
    cropy= int(x*perc)
    startx = x//2 - cropx//2
    starty = y//2 - cropy//2    
    return matrix[starty:starty+cropy, startx:startx+cropx]

def crop_array(array,perc=0.3):
    out=[]
    for arr in array:
        x = arr.shape[0]
        cropx= int(x*perc)
        startx = x//2 - cropx//2
        arr = arr[startx:startx+cropx]
        out.append(arr)
    return out


percent=0.3

mean_matrix_r = crop_center(mean_matrix,percent )
variance_matrix_r = crop_center(variance_matrix, percent)
xy_r=crop_array(xy)

start_r =-1.0  + (1-percent)
leng_r = 2.0 - 2.0 * (1-percent)

# %% Centered


fig = show_mean_varaince_LOG_2d_3d(mean_matrix_r, variance_matrix_r, xy_r, BW= False, start=start_r, leng=leng_r)

formula=setup_dict['eq_formula']
specifc=setup_dict['eq_formula2']
title_txt = f'Center solutions Over:{N} \n {formula}  \n {specifc}'  #
fig.suptitle(title_txt)
fig.tight_layout()
fig.show()


file_name = f'all_Log_solutionsMeanVar_center_'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh']
fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")


# %% Centered





from cdf_project.data_exp.explore_utils_fileds import show_mean_varaince_2d_3d
fig = show_mean_varaince_2d_3d(mean_matrix_r, variance_matrix_r, xy_r, BW= False, start=start_r, leng=leng_r)


formula=setup_dict['eq_formula']
specifc=setup_dict['eq_formula2']
title_txt = f'Center solutions Over:{N} \n {formula}  \n {specifc}'  #
fig.suptitle(title_txt)
fig.tight_layout()
fig.show()


file_name = f'all_solutionsMeanVar_center_'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh']
fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")



fig = show_mean_varaince_2d_3d(mean_matrix_r, variance_matrix_r, xy_r, BW= True, start=start_r, leng=leng_r)


formula=setup_dict['eq_formula']
specifc=setup_dict['eq_formula2']
title_txt = f'Center solutions Over:{N} \n {formula}  \n {specifc}'  #
fig.suptitle(title_txt)
fig.tight_layout()
fig.show()


file_name = f'all_solutionsMeanVar_center_'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + '_bw'
fig.savefig(report_folder/(file_name+'.jpg'), bbox_inches="tight")




# %%  Salve overall intresting points....

from cdf_project.data_etl.etl_utils import array_from_mesh
xy=array_from_mesh(setup_dict['dimensions'], input_mesh)


points_dict={}
points_dict['setup_dict']=setup_dict
points_dict['center_values_array']=center_values
points_dict['seeds_array'] = seeds_array
points_dict['center_value_x']= f"x({xy[0][nxy//2]:.4f} @ {nxy//2})"
points_dict['center_value_y']= f"y({xy[1][nxy//2]:.4f} @ {nxy//2})"

#%% Save or create
sol_hash_name = setup_dict['sol_hash_name']
point_dir = f'points_{sol_hash_name}'
(point_data_dir / point_dir).mkdir(exist_ok=True) #parents=True, 

point_file_name=f'points_{sol_hash_name}.plk'
points_file = point_data_dir / point_dir / point_file_name

#%%  if points exists read evals
if points_file.is_file():
    print('Exist, loading')
    with open(points_file, 'rb') as f:
        solutions_dict = pickle.load(f)
else:
    print('Not exist, saving')
    with open(points_file, 'wb') as f:
        pickle.dump(points_dict, f)

# %%
