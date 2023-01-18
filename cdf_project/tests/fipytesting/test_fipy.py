
#%% 1.IMprts 
from cdf_project.setup_testSet import setup_dict

from cdf_project.data_etl.etl_utils import gen_fields, gen_uniform_delta


setup_dict['number_of_cells'] = 1000
nxy=setup_dict['number_of_cells']
dim = 2
#0.9 0.25  15
corel=0.5
variance= 0.01
setup_filed={'dimensions' : dim,
'correlation_coeff': corel,
'variance':variance}

dict_out = \
    gen_fields(dims= dim, number_of_cells=nxy, input_dict=setup_filed,\
         ens_num=4, enstart= 0, log_normal=False, \
        non_uniform='uniform', interval_len = 2.0, interval_start= -1.0 )

fields = dict_out['fields_array']
mesh =  dict_out['mesh']


delta_space = gen_uniform_delta(nxy)  
setup_dict['grid_spacing'] = delta_space

#%% 1.IMprts 
from cdf_project.feature_eng.solver4stationary_util import setup_phi_0, get_dirac, get_equation, swipe_solve
variables_dict = setup_dict


field_number= 1
flatten_field = fields[field_number].flatten()






#%% 2-3.  Set initial conditions and boundary conditions
phi = setup_phi_0(variables_dict, mesh)

#%% 4-5.  set Param AND equation
#equ = get_equation(variables_dict, mesh, flatten_field, phi)
import fipy as fi
fpnp = fi.numerix
    
#source = 'eFdelta'; D = 1.0; convect = True 
source = 'delta'; D = 0.0; convect = False 


F= fi.CellVariable(mesh=mesh, rank=0)
if source == 'delta':
    F[:] =  get_dirac(0.0, variables_dict, mesh).value #- (1-x) 
    sourceCoeff = F
elif source == 'Fdelta':
    F[:] =  flatten_field #- (1-x) 
    sourceCoeff = F  * get_dirac(0.0, variables_dict, mesh).value
elif source == 'eFdelta':
    F[:] =  flatten_field #- (1-x) 
    sourceCoeff = fpnp.exp(-F) * get_dirac(0.0, variables_dict, mesh).value
else:
    raise NameError('Not a valid source type')


if D == 0.0:
    D=fi.CellVariable(mesh=mesh, rank=0, value=fpnp.exp(flatten_field))
    #D[:] =  fpnp.exp(flatten_field)


U= fi.CellVariable(mesh=mesh, rank=0)
U[:] =  flatten_field #- (1-x)
#TODO gradient of U
grad_U =  U.grad  #faceGrad  grad
if convect:
    equ =  (fi.DiffusionTerm(coeff=D, var=phi)
            + fi.ExponentialConvectionTerm(coeff=grad_U, var=phi)
            - sourceCoeff)
elif not convect:
    equ =  (fi.DiffusionTerm(coeff=D, var=phi)
            - sourceCoeff)   # == 0.0
else:
    raise NameError(f'Not contemplated setup convect:{convect}')


#%% 7. Solving PDE
dim = variables_dict['dimensions']
phi_tt, _  =swipe_solve(phi, equ, dim,  verbose=False)










# %%
from cdf_project.data_exp.explore_utils_sols import show_sinlge_2d_sol_filed
dimens = setup_dict['dimensions']
nxy = setup_dict['number_of_cells']

values_matrix = phi_tt.reshape((nxy, nxy))
field=fields[0]
import numpy as np
ind = [nxy//2,nxy//2]


from cdf_project.data_etl.etl_utils import array_from_mesh
xy=array_from_mesh(2, mesh)

fig, ax = show_sinlge_2d_sol_filed(field, values_matrix, xy, ind)

# %%
import os
solver =  os.environ.get('FIPY_SOLVERS')


title_txt = f'fld{field_number} GMRE delta:{source} convect:{convect} cor{corel} variance {variance} mesh{nxy}'  #
fig.suptitle(title_txt)
fig.tight_layout()

from cdf_project.config import reportings_dir


file_name = f'cor{int(corel*100)} variance {int(variance*100)} fld{field_number}  mesh{nxy} _   delta:{source} convect:{convect} GMRE{solver}'  #
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
