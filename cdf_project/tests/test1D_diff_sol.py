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

#%%  Import pakages
from cdf_project.config import data_dir, reportings_dir
from cdf_project.setup_testSet import setup_dict
from cdf_project.tests.tests_utils import plot_conv_acc, calulate_all_errors, plot_final_funcs_x


import numpy as np
import matplotlib.pyplot as plt
import fipy as fp

##%%
#import matplotlib.pyplot as plt
#plt.style.use('seaborn-notebook')








#%%  Just debug test
from cdf_project.data_etl.etl_utils import gen_nonuniform_segment_delta, gen_nonuniform_chebyshev_delta, gen_x_nnunif_mesh, gen_x_mesh, gen_uniform_delta
from cdf_project.feature_eng.solver_util import de_f_time

number_of_cells=1200
delta_space = gen_uniform_delta(number_of_cells)
mesh, x =gen_x_mesh(number_of_cells, delta_space)
setup_dict['grid_spacing'] = delta_space
setup_dict['center_spacing'] = delta_space


#delta_vect, delta_space = gen_nonuniform_chebyshev_delta(number_of_cells)
#mesh, x =gen_x_nnunif_mesh(delta_vect)
#setup_dict['grid_spacing'] = max(delta_vect)
#setup_dict['center_spacing'] = max(delta_vect)


timeDuration = 0.04
three_time_phi = de_f_time(1, setup_dict, mesh, None, timeDuration)  #False  True
assert(len(three_time_phi[1][-1])==number_of_cells)




#%%  Read fileds and solution
#Only generate single - determinisctic solution solution
timeDuration = 0.04


#%%  save a VIDEO!


#%%  Error testing 
def find_error(mesh, numerc_sol, t):
        global setup_dict
        D= setup_dict['D']
        fpnp = fp.numerix
        analyt_sol = fp.CellVariable(name="analytical solution", mesh=mesh, value=0.)
        den = 4*fpnp.pi*D*t #TODO   4*fpnp.pi*t o D?
        X = mesh.x #mesh.cellCenters[0].value
        analyt_sol[:] = (1/ (den**0.5)) * fpnp.exp(-(X-0.5)**2/(4*D*t))
        error = (fpnp.sum((numerc_sol - analyt_sol)**2)/len(numerc_sol))**(0.5)
        return error.value, analyt_sol.value





#%%  Plot single errors
fig =  plot_final_funcs_x(150*5, find_error, timeDuration)




#%%  Generate num for convergence accuracy
nxvec = np.arange(720, 750)
errors_u, errors_nu_seg, errors_nu_chb = calulate_all_errors(nxvec, find_error, timeDuration)

#%%  Plot
fig = plot_conv_acc(timeDuration, nxvec[15:30], errors_u[15:30], errors_nu_seg[15:30], errors_nu_chb[15:30])






















# %%
