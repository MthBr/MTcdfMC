

#%%
from cdf_project.data_etl.etl_utils import gen_fields

nxy=8
dim = 2
dict_out = \
    gen_fields(dims= dim, number_of_cells=nxy, correl=0.9,\
         ens_num=4, enstart= 0, log_normal=False, \
        non_uniform='uniform', interval_len = 2.0, interval_start= -1.0 )

fields = dict_out['fields_array']
mesh =  dict_out['mesh']


#%%
from cdf_project.data_exp.explore_utils_fileds import plot4fileds
fig, ax = plot4fileds(fields, dim)




# %%
from cdf_project.data_exp.explore_utils_fileds import filed2d_in_3d
from cdf_project.data_etl.etl_utils import array_from_mesh
d_vector_mesh = array_from_mesh(dim, mesh)
fig, ax = filed2d_in_3d(fields[0], d_vector_mesh)
# %%

import fipy as fi

U= fi.CellVariable(mesh=mesh, rank=1)
flatten_field =   fields[1].flatten()
U[:] =  flatten_field #- (1-x) 
grad_U =  U.grad #faceGrad  grad

matrix = [grad_U[0][0].value.reshape((nxy, nxy)),\
    grad_U[0][1].value.reshape((nxy, nxy)),\
      grad_U[1][0].value.reshape((nxy, nxy)),\
          grad_U[1][1].value.reshape((nxy, nxy))]  



fig, ax = filed2d_in_3d(matrix[1], d_vector_mesh)


fig, ax = plot4fileds(matrix, dim)



# %%
import fipy as fi
from cdf_project.feature_eng.delta_func_utils import delta_func
from cdf_project.data_etl.etl_utils import gen_uniform_delta
flatten_field =   fields[1].flatten()
F= fi.CellVariable(mesh=mesh, rank=0)
fpnp = fi.numerix
F[:] =  flatten_field #- (1-x) 

delta_space  = gen_uniform_delta(number_of_cells=nxy, length=2.0)
dirac = fi.CellVariable(name="dirac", mesh=mesh,  value=0.)
appx_type = '2*-cos'
center = 0.0
dirac[:] = delta_func(mesh, center, h=1.4*delta_space,\
            dim=dim, approx_type= appx_type)


sourceCoeff = fpnp.exp(-F) * dirac.value



# %%
matrixF = sourceCoeff.value.reshape((nxy, nxy))


fig, ax = filed2d_in_3d(matrixF, d_vector_mesh)


fig, ax = plot4fileds([matrixF], dim)
# %%

# %%
