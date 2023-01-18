#%% Clean Memory and import pakages
for i in list(globals().keys()):
    if(i[0] != '_'):
        exec('del {}'.format(i))

from cdf_project.config import data_dir
from cdf_project.setup_testSet import setup_dict

import pickle

#%%  Dataset paths
raw_data_dir = data_dir / 'raw'
point_data_dir = data_dir / 'models'


#%% Save or create
sol_hash_name = setup_dict['sol_hash_name']
point_dir = f'points_{sol_hash_name}'
point_file_name=f'points_{sol_hash_name}.plk'
points_file = point_data_dir / point_dir / point_file_name

#%%  if points exists read evals
if points_file.is_file():
    print('Exist, loading')
    with open(points_file, 'rb') as f:
        points_dict = pickle.load(f)
else:
    print(f'Not exist!!!!{points_file}')




#%%

setup_dict= points_dict['setup_dict']
center_values = points_dict['center_values_array']
seeds_array = points_dict['seeds_array']
x_label = points_dict['center_value_x']
y_label = points_dict['center_value_y']


#%% 7. CDF plot

from util_stats_pdf import cdf_1

fig = cdf_1(center_values, x_label, y_label)

fig.suptitle(f"{setup_dict['eq_formula']} \n {setup_dict['eq_formula2']}")

fig.tight_layout()
fig.show()

from cdf_project.config import reportings_dir
file_name = f'CDF1'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq']
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



#%% 7. CDF plot
from util_stats_pdf import cdf_2

fig = cdf_2(center_values, x_label, y_label)
fig.suptitle(f"{setup_dict['eq_formula']} \n {setup_dict['eq_formula2']}")

fig.show()
fig.tight_layout()
from cdf_project.config import reportings_dir
file_name = f'CDF2'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq']
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



#%% 7. cumulative:  1.mean - 2.std - 3.   plot
#https://stackoverflow.com/questions/53005146/calculating-the-cumulative-mean-in-python

probe = center_values # almeno 2

mc_runs = len(probe)+1


from util_stats_pdf import cumulative

fig = cumulative(probe, mc_runs)

fig.suptitle(f"{setup_dict['eq_formula']} \n {setup_dict['eq_formula2']}")

fig.tight_layout()
fig.show()
#plt.savefig(f'moments_{cr_txt}.jpg', dpi=300)


fig.show()
fig.tight_layout()
from cdf_project.config import reportings_dir
file_name = f'order'+ setup_dict['compactxt_field'] + setup_dict['compactxt_mesh'] + setup_dict['compactxt_eq']
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

# %%
