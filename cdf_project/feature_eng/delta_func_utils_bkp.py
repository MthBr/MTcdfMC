# -*- coding: utf-8 -*-
"""
Version 1, smoothing the Dirac delta function 
 as is done with diffusion interface methods 
 see equations 11, 12 and 13:
 https://doi.org/10.1016/j.jcp.2004.09.018
AND
Eq 13-17 in


@author: enzo
"""

#%% import pakages
from logging import error
import fipy as fi
import numpy as np


from cdf_project.custom_funcs import  get_logger, benchmark
from cdf_project.config import log_dir
global logger
if __name__ == '__main__':
    logger = get_logger(log_dir/'delta_functs.log', 'general')
else:
    logger = get_logger()




#%% Load functions
global max_delta
max_delta = 150.1449498535622
global space_offset
space_offset=10**-7

@benchmark
def delta_old(mesh, delta_space, dim):
    #Initial Impulse condition
    phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    delta= delta_space - delta_space*space_offset #10**-3     #delta_space - delta_space*0.3
    where_delta= (abs(mesh.cellCenters[0]-(0.5+space_offset)) < delta)
    if dim >= 2 : where_delta = where_delta & (abs(mesh.cellCenters[1]-0.5) < delta)
    if dim >= 3 : where_delta = where_delta & (abs(mesh.cellCenters[2]-0.5) < delta)
    phi.setValue(1.*max_delta, where=where_delta )  #dxy*nxy//2
    assert np.sum(np.array(phi.value) > 0) == 2**dim
    return phi

@benchmark
def delta_c_func_my(mesh, x0, epsilon, dim):
    # this implementation is FASTER!???
    fpnp = fi.numerix
    x = mesh.cellCenters[0]-x0 # mesh.x-x0
    where_delta= (abs(x) < epsilon)
    delt_values = (1 + fpnp.cos(fpnp.pi * x / epsilon)) / 2 / epsilon
    if dim >= 2 : 
        y = mesh.cellCenters[1]-x0
        where_delta = where_delta & (abs(y) < epsilon)
        delt_values = delt_values *  (1. + fpnp.cos(fpnp.pi * y / epsilon)) / 2 / epsilon
    if dim >= 3 : 
        z = mesh.cellCenters[2]-x0
        where_delta = where_delta & (abs(z) < epsilon)
        delt_values = delt_values *  (1. + fpnp.cos(fpnp.pi * z / epsilon)) / 2 / epsilon
    delt_values = where_delta * delt_values
    logger.debug(f'MY test Delta cos Dirac non negative: {np.sum(np.array(delt_values.value) > 0)}')
    return delt_values.value

@benchmark
def delta_c_func(mesh, x0, epsilon, dim):
    import numpy as np
    fpnp = fi.numerix

    x = mesh.cellCenters[0]-x0  # mesh.x- x0
    where_delta_x= (abs(x) < epsilon)
    delt_values = where_delta_x * \
        (1. + fpnp.cos(fpnp.pi * x / epsilon)) / 2 / epsilon

    if dim >= 2 : 
        y = mesh.y-x0  #mesh.cellCenters[1]-x0
        where_delta_y = (abs(y) < epsilon)
        delt_y =  where_delta_y * (1. + fpnp.cos(fpnp.pi * y / epsilon)) / 2 / epsilon
        delt_values = delt_values *  delt_y
    if dim >= 3 : 
        z = mesh.z-x0  #mesh.cellCenters[2]-x0
        where_delta_z = (abs(z) < epsilon)
        delt_z =  where_delta_z * (1. + fpnp.cos(fpnp.pi * z / epsilon)) / 2 / epsilon
        delt_values = delt_values *  delt_z

    logger.debug(f'Delta cos Dirac non negative: {np.sum(np.array(delt_values.value) > 0)}')
    return delt_values.value


def delta_cos_func(mesh, x0, h, dim, m=2):
    #4-point cosine function
    return delta_func(mesh, x0, h, dim, approx_type= '2-cos' ,m=2)

def delta_l_func(mesh, x0, h, dim, m=1):
    #2-point hat function
    return delta_func(mesh, x0, h, dim, approx_type= '1-l*' ,m=1)



def delta_cubic_func(mesh, x0, h, dim, m=2):
    return delta_func(mesh, x0, h, dim, approx_type= 'Cubic' ,m=2)


def delta_LL_func(mesh, x0, h, dim, m=2):
    return delta_func(mesh, x0, h, dim, approx_type= 'LL' ,m=2)

# Notes
# cos : 15 punt
# l : 3 - 5 -7


#%% Define known deltas





#%%


def approx(approx_type, x, h, m):
    fpnp = fi.numerix

    if approx_type == '1-l':
        #2-point hat function
        moment = m = 1
        epsilon = m*h
        where_delta= (abs(x) <= epsilon)
        phi = 1./(m**2) * np.minimum(x/h + m, m - x/h)
    elif approx_type == '2-l':
        #4-point hat function
        moment = m = 2
        epsilon = m*h
        where_delta= (abs(x) <= epsilon)
        phi = 1./(m**2) * np.minimum(x/h + m, m - x/h)
    elif approx_type == '1*-l':
        #smoothed 2-point hat function
        moment = m = 1
        epsilon = (m+0.5)*h
        xi = abs(x/h).value
        where_delta= (abs(xi) <= epsilon/h) #TODO does not work!!!
        phi =np.where(xi <= 0.5, \
            0.75 - xi**2, \
            1.125 - 1.5*xi + 0.5* xi**2)
    elif approx_type == '1-l*':
        #smoothed 2-point hat function
        moment = m = 1
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)
        xi = abs(x/h).value
        phi =np.where(xi <= 0.5, \
            0.75 - xi**2, \
            1.125 - 1.5*xi + 0.5* xi**2)
    elif approx_type == '1-cos':
        #4-point cosine function
        moment = m = 1
        epsilon = m*h
        where_delta= (abs(x) <= epsilon)
        ra = abs(x/h).value
        phi = 1./2 * (1+ fpnp.cos(fpnp.pi * ra))
    elif approx_type == '2-cos':
        #4-point cosine function
        moment = m = 2
        epsilon = m*h
        where_delta= (abs(x) <= epsilon)
        phi = 1./4 * (1+ fpnp.cos(fpnp.pi * x / (2*h)))
    elif approx_type == 'cos_altern':
        moment = m = 2
        epsilon = m*h
        where_delta= (abs(x) < epsilon)
        phi = (1. + fpnp.cos(fpnp.pi * x / epsilon)) / 2 / epsilon
    elif approx_type == '2*-cos':
        #smoothed 4-point cosine function - Class C^2
        moment = m = 2
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)
        ra = abs(x/h).value
        r = (x/h).value
        Pi=fpnp.pi
        Pi4 = fpnp.pi/4
        Pi2 = fpnp.pi/2
        phi1 =np.where(ra <= 1.5, \
            1/4 + 1/(2*Pi)* fpnp.sin(Pi2*r+Pi4) - 1/(2*Pi)*fpnp.sin(Pi2*r-Pi4)  , \
                0.0 )
        phi2 =np.where((1.5 < ra) &  (ra<= 2.5), \
            5/8 - ra/4 - 1/(2*Pi) *fpnp.sin(Pi2*ra-Pi4) , \
             0.0)
        phi = phi1+phi2
    elif approx_type == '2*-cosTEST3':
        #smoothed 4-point cosine function - Class C^2
        moment = m = 2
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)
        ra = abs(x/h).value
        r = (x/h).value
        Pi=fpnp.pi
        Pi4 = fpnp.pi/4
        Pi2 = fpnp.pi/2

        mask1 = (ra <= 1.5)
        phi1 = np.zeros_like(ra)
        r_m = r[mask1]
        phi1[mask1] =  1/4 + 1/(2*Pi)* fpnp.sin(Pi2*r_m+Pi4) - 1/(2*Pi)*fpnp.sin(Pi2*r_m-Pi4)
        #phi1[~mask1] = 0.0

        mask2 = ((1.5 < ra) &  (ra<= 2.5))
        phi2 = np.zeros_like(ra)
        ra_m = ra[mask2]
        phi2[mask2] = 5/8 - ra_m/4 - 1/(2*Pi) *fpnp.sin(Pi2*ra_m-Pi4)
        #phi2[~mask2] = 0.0
        phi = phi1+phi2
    elif approx_type == '2*-cosTEST':
        #smoothed 4-point cosine function - Class C^2
        moment = m = 2
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)
        ra = abs(x/h).value
        r = (x/h).value
        Pi=fpnp.pi
        phi =np.where(ra <= 1.5, \
            1/(4*Pi)*( Pi + 2* fpnp.sin(Pi/4*(2*r+1)) - 2*fpnp.sin(Pi/4*(2*r-1)) )   , \
            -1/(8*Pi)*(-5*Pi+2*Pi*ra  + 4*fpnp.sin(Pi/4*(2*ra-1))) )
    elif approx_type == '2*-cosTEST2':
        #smoothed 4-point cosine function - Class C^2
        moment = m = 2
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)
        ra = abs(x/h).value
        r = (x/h).value
        Pi=fpnp.pi
        Pi4 = fpnp.pi/4
        phi1 =np.where(ra <= 1.5, \
            1/(4*Pi)*( Pi + 2* fpnp.sin(Pi4*(2*r+1)) - 2*fpnp.sin(Pi4*(2*r-1)) ) , \
                0.0 )
        phi2 =np.where((1.5 < ra) &  (ra<= 2.5), \
            -1/(8*Pi)*(-5*Pi+2*Pi*ra  + 4*fpnp.sin(Pi4*(2*ra-1))) , \
             0.0)
        phi = phi1+phi2
    elif approx_type == '3*f-arc':
        #smoothed 3-point function - Class C^2
        moment = m = 3
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)
        ra = abs(x/h).value
        sq3 = fpnp.sqrt(3)
        Pi=fpnp.pi
        sq3P108=sq3*Pi/108
        phi1 =np.where(ra <= 1.0, \
            17/48 + sq3P108 + ra/4 - ra**2/4 + (1-2*ra)/16*(fpnp.sqrt(-12*ra**2+12*ra+1)) - sq3/12*fpnp.arcsin(sq3/2*(2*ra-1)) , \
                0.0 )
        phi2 =np.where((1.0 < ra) &  (ra<= 2), \
            55/48 - sq3P108 -13*ra/12 +ra**2/4 + (2*ra-3)/48*(fpnp.sqrt(-12*ra**2+36*ra-23)) + sq3/36*fpnp.arcsin(sq3/2*(2*ra-3))      , \
             0.0)
        phi = phi1+phi2
    elif approx_type == '4*f':
        #smoothed 4-point piecewise function - Class C^2
        moment = m = 4
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)
        ra = abs(x/h).value
        sq2 = fpnp.sqrt(2)
        Pi=fpnp.pi
        phi1 =np.where(ra <= 0.5, \
            3/8 + Pi/32 - ra**2/4 , 0.0 )
        phi2 =np.where((0.5 < ra) &  (ra<= 1.5), \
            1/4+ (1-ra)/8*(fpnp.sqrt(-2 + 8*ra -4**2*ra)) -1/8*fpnp.arcsin(sq2*(ra-1)), \
             0.0)
        phi3 =np.where((1.5 < ra) &  (ra<= 2.5), \
            17/16-Pi/64-3*ra/4+ra**2/8+(ra-2)/16*(fpnp.sqrt(-14+16*ra-4**2*ra)) + 1/16*fpnp.arcsin(sq2*(ra-2)), \
             0.0)        
        phi = phi1+phi2+phi3
    elif approx_type == '4f':
        #smoothed 4-point piecewise function - Class C^2
        moment = m = 2
        epsilon = (m+0.5)*h
        where_delta= (abs(x) <= epsilon)

        ra = abs(x/h).value
        sq2 = fpnp.sqrt(2)
        Pi=fpnp.pi

        mask1 = (ra <= 0.5)
        phi1 = np.zeros_like(ra)
        ra_m = ra[mask1]
        phi1[mask1] = 3/8 + Pi/32 - ra_m**2/4
        phi1[~mask1] = 0.0

        mask2 = ((0.5 < ra) &  (ra<= 1.5))
        phi2 = np.zeros_like(ra)
        ra_m = ra[mask2]
        phi2[mask2] = 1/4+ (1-ra_m)/8*(np.sqrt(-2 + 8*ra_m -4*ra_m**2)) -1/8*np.arcsin(sq2*(ra_m-1))
        phi2[~mask2] = 0.0

        mask3 = ((1.5 < ra) &  (ra<= 2.5))
        phi3 = np.zeros_like(ra)
        ra_m = ra[mask3]
        phi3[mask3] = 17/16-Pi/64-3*ra_m/4+ra_m**2/8+(ra_m-2)/16*(np.sqrt(-14+16*ra_m-4*ra_m**2)) + 1/16*np.arcsin(sq2*(ra_m-2))
        phi3[~mask3] = 0.0 
        phi = phi1+phi2+phi3
    elif approx_type == '2-Cubic':
        moment = m = 2
        epsilon = m*h
        where_delta= (abs(x) <= epsilon)
        xi = abs(x/h).value
        phi =np.where(xi <= 1, \
            1. - 0.5*xi - xi**2 + 0.5 * xi**3 , \
            1. - 11./6.*xi + xi**2 - 1./6.* xi**3)
    elif approx_type == '2-LL':
        moment = m = 2
        epsilon = m*h
        where_delta = (abs(x) <= epsilon)
        xi = abs(x/h).value
        phi =np.where(xi <= 1, \
            1./12. * (14-15*xi), \
            1./12. * (2-xi) )
    else:
        raise NameError(f'Not a valid approx_type:  {approx_type}')
    
      
    delta =  where_delta * 1./h * phi
    # h == delta_space
    if h*moment < 9.9*10**-3:
        logger.debug(f'({moment}=m)*epsilon= {h*moment} < {9.9*10**-3}')
    else:
        logger.critical(f'({moment}=m)*epsilon= {h*moment} > {9.9*10**-3}')
    #assert(h*moment < 9.9*10**-3)

    return delta

def delta_func(mesh, x0, h, dim, approx_type, m=2):
    import numpy as np
    x = mesh.cellCenters[0]-x0  # mesh.x- x0
    delt_values = approx(approx_type, x, h, m)
    logger.debug(f'did 1D layer')
    if dim >= 2 : 
        logger.debug(f'2D layer')
        y = mesh.y-x0  #mesh.cellCenters[1]-x0
        delt_values = delt_values * approx(approx_type, y, h, m)
    if dim >= 3 : 
        logger.debug(f'add 3D layer')
        z = mesh.z-x0  #mesh.cellCenters[2]-x0
        delt_values = delt_values * approx(approx_type, z, h, m)
    logger.debug(f'Delta {approx_type} Dirac non negative: {np.sum(np.array(delt_values.value) > 0)}')
    return delt_values.value


def delta_radial_func(mesh, x0, h, dim, approx_type, m=2):
    #Initial Impulse condition
    # delta_space = h
    xyz= abs(mesh.cellCenters[0]-x0)
    if dim >= 2 : 
        xyz = np.sqrt((mesh.cellCenters[0]-x0)**2 + (mesh.cellCenters[1]-x0)**2)
    if dim >= 3 : 
        xyz = np.sqrt((mesh.cellCenters[0]-x0)**2 + (mesh.cellCenters[1]-x0)**2 + (mesh.cellCenters[2]-x0)**2)

    delt_values =  approx(approx_type, xyz, h, m)

    logger.debug(f'Delta {approx_type} Dirac non negative: {np.sum(np.array(delt_values.value) > 0)}')
    return delt_values.value











@benchmark
def test_delta_func(mesh, x0, h, dim, approx_type, m=2):
    for _ in range(3):
        t = delta_func(mesh, x0, h, dim, approx_type, m)
    return t






#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import array_from_mesh, gen_xyz_nnunif_mesh, gen_nonuniform_segment_delta
    dim = 3
    nxy = number_of_cells = 91
    delta_vect, delta_space = gen_nonuniform_segment_delta(number_of_cells, length=2)
    mesh =gen_xyz_nnunif_mesh(delta_vect, start=-1)
    xyz=array_from_mesh(3, mesh, len(delta_vect))
    epsilon = delta_space *1.4
    TEST =delta_func(mesh, 0.0, epsilon, dim, '2*-cosTEST') #TODO Just for testing Purpose
    values = delta_func(mesh, 0.0, epsilon, dim, '2*-cosTEST2')
    np.testing.assert_array_equal(TEST, values)
    TEST_3d=TEST.reshape((nxy, nxy, nxy))
    values_3d = values.reshape((nxy, nxy, nxy))

    ind = np.unravel_index(np.argmax(values_3d, axis=None), values_3d.shape)
    print(ind)
    print([nxy//2,nxy//2,nxy//2])

    from cdf_project.data_exp.explore_utils_sols import compare_show_3d
    fig, axs = compare_show_3d(values_3d, TEST_3d, xyz, loc = ind, uniform= False)
    fig.show()






#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import array_from_mesh, gen_xyz_mesh, gen_uniform_delta
    dim = 3
    nxy = number_of_cells = 125
    delta_space = gen_uniform_delta(number_of_cells, length=2)
    mesh =gen_xyz_mesh(number_of_cells, delta_space, start=-1)
    xyz=array_from_mesh(3, mesh)
    epsilon = delta_space *1.4
    TEST =delta_func(mesh, 0.0, epsilon, dim, '2*-cos') #TODO Just for testing Purpose
    values = delta_func(mesh, 0.0, epsilon, dim, '2*-cosTEST2')
    np.testing.assert_array_equal(TEST, values)
    TEST_3d=TEST.reshape((nxy, nxy, nxy))
    values_3d = values.reshape((nxy, nxy, nxy))

    from cdf_project.data_exp.explore_utils_sols import compare_show_3d
    fig, axs = compare_show_3d(values_3d, TEST_3d, xyz, loc = [nxy//2,nxy//2,nxy//2])
    fig.show()



#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_xyz_mesh, gen_uniform_delta
    dim = 3
    nxy = number_of_cells = 150
    delta_space = gen_uniform_delta(number_of_cells, length=2)
    mesh =gen_xyz_mesh(number_of_cells, delta_space, start=-1)
    epsilon = delta_space *1.4

    test_delta_func(mesh, 0.0, epsilon, dim, '2*-cos') #TODO Just for testing Purpose
    
    test_delta_func(mesh, 0.0, epsilon, dim, '2*-cosTEST2')





#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import array_from_mesh, gen_xyz_mesh, gen_uniform_delta
    dim = 3
    nxy = number_of_cells = 125
    delta_space = gen_uniform_delta(number_of_cells, length=2)
    mesh =gen_xyz_mesh(number_of_cells, delta_space, start=-1)
    xyz=array_from_mesh(3, mesh)
    epsilon = delta_space *1.4
    TEST =delta_func(mesh, 0.0, epsilon, dim, '3*f-arc') #TODO Just for testing Purpose
    values = delta_func(mesh, 0.0, epsilon, dim, '2-Cubic') # 2*-cos 2-Cubic
    #np.testing.assert_array_equal(TEST, values)
    TEST_3d=TEST.reshape((nxy, nxy, nxy))
    values_3d = values.reshape((nxy, nxy, nxy))

    from cdf_project.data_exp.explore_utils_sols import compare_show_3d
    fig, axs = compare_show_3d(values_3d, TEST_3d, xyz, loc = [nxy//2,nxy//2,nxy//2])
    fig.show()






#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_x_mesh, gen_uniform_delta
    dim = 1
    number_of_cells = 300
    delta_space = gen_uniform_delta(number_of_cells, 2)
    mesh, x =gen_x_mesh(number_of_cells, delta_space, -1)
    epsilon = delta_space *1.4
    TEST =delta_func(mesh, 0.0, epsilon, dim, '2*-cosTEST') #TODO Just for testing Purpose
    values = delta_func(mesh, 0.0, epsilon, dim, '2*-cosTEST2')
    np.testing.assert_array_equal(TEST, values)

    import matplotlib.pyplot as plt
    plt.plot(x[0], TEST, label="1")
    plt.plot(x[0], values, linestyle="--", label="2")
    plt.show()















#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_x_mesh, gen_uniform_delta
    dim = 1
    number_of_cells = 300
    delta_space = gen_uniform_delta(number_of_cells, 2)
    mesh, x =gen_x_mesh(number_of_cells, delta_space, -1)
    epsilon = delta_space *1.4
    TEST =delta_func(mesh, 0.0, epsilon, dim, '1-l*') #TODO Just for testing Purpose
    values = delta_func(mesh, 0.0, epsilon, dim, '1*-l')
    np.testing.assert_array_equal(TEST, values)
    import matplotlib.pyplot as plt
    plt.plot(x[0], TEST)
    plt.plot(x[0], values)
    plt.show()


#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_x_mesh, gen_uniform_delta
    dim = 1
    number_of_cells = 150
    delta_space = gen_uniform_delta(number_of_cells)
    mesh, _ =gen_x_mesh(number_of_cells, delta_space)
    mesh = (mesh - [[0.5]])*4
    x = [mesh.x.value]
    epsilon = delta_space *number_of_cells

    phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    #phi_t = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    #phi_t[:] =delta_c_func_my(mesh, 0, epsilon*2, dim) #TODO Just for testing Purpose
    #phi[:] = delta_LL_func(mesh, 0, epsilon, dim) # delta_cos_func delta_l_func
    phi[:] = delta_cos_func(mesh, 0, epsilon, dim) # delta_cubic_func delta_LL_func

    #np.testing.assert_array_equal(phi_t.value, phi.value)
    import matplotlib.pyplot as plt
    #plt.plot(x[0], phi_t.value)
    plt.plot(x[0], phi.value)
    plt.show()
    print(min(phi.value), max(phi.value))


#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_x_mesh, gen_uniform_delta
    dim = 1
    number_of_cells = 400
    delta_space = gen_uniform_delta(number_of_cells)
    mesh =gen_x_mesh(number_of_cells, delta_space)
    x = array_from_mesh(1, mesh)
    epsilon = delta_space *5

    phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    #phi_t = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    #phi_t[:] =delta_c_func_my(mesh, 0, epsilon*2, dim) #TODO Just for testing Purpose
    phi[:] = delta_cos_func(mesh, 0.5, epsilon, dim) # delta_cos_func delta_l_func
    #phi[:] = delta_cubic_func(mesh, 0.5, epsilon, dim) # delta_cubic_func delta_LL_func

    #np.testing.assert_array_equal(phi_t.value, phi.value)
    import matplotlib.pyplot as plt
    #plt.plot(x[0], phi_t.value)
    plt.plot(x[0], phi.value)
    plt.show()
    print(min(phi.value), max(phi.value))

#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_x_mesh, gen_uniform_delta
    dim = 1
    number_of_cells = 125
    delta_space = gen_uniform_delta(number_of_cells)
    mesh =gen_x_mesh(number_of_cells, delta_space)
    x = array_from_mesh(1, mesh)
    epsilon = delta_space *10.
    TEST =delta_c_func_my(mesh, 0.5, epsilon, dim) #TODO Just for testing Purpose
    values = delta_c_func(mesh, 0.5, epsilon, dim)
    np.testing.assert_array_equal(TEST, values)
    import matplotlib.pyplot as plt
    plt.plot(x[0], TEST)
    plt.plot(x[0], values)
    plt.show()




#%%
if __name__ == '__main__':
    dim = 1
    number_of_cells = 5
    delta_space = 1
    mesh =gen_x_mesh(number_of_cells, delta_space)
    mesh = mesh - [[2.5]]
    x = [mesh.x.value]
    epsilon = delta_space *3.

    phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    phi_t = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    phi_t[:] =delta_c_func_my(mesh, 0, epsilon, dim) #TODO Just for testing Purpose
    phi[:] = delta_c_func(mesh, 0, epsilon, dim)

    np.testing.assert_array_equal(phi_t.value, phi.value)
    import matplotlib.pyplot as plt
    plt.plot(x[0], phi_t.value)
    plt.plot(x[0], phi.value)
    plt.show()








#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_xy_mesh, gen_uniform_delta
    dim = 2
    number_of_cells = 125
    delta_space = gen_uniform_delta(number_of_cells)
    mesh =gen_xy_mesh(number_of_cells, delta_space)
    epsilon = delta_space *2.
    TEST =delta_c_func_my(mesh, 0.5, epsilon, dim) #TODO Just for testing Purpose
    values = delta_c_func(mesh, 0.5, epsilon, dim)
    np.testing.assert_array_equal(TEST, values)
    import matplotlib.pyplot as plt
    plt.plot(x[0], TEST)
    plt.plot(x[0], values)
    plt.show()

#%%






if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_xyz_mesh, gen_uniform_delta
    dim = 3
    nxy = number_of_cells = 125
    delta_space = gen_uniform_delta(number_of_cells, length=2)
    mesh =gen_xyz_mesh(number_of_cells, delta_space, start=-1)
    
    epsilon = delta_space *1.4
    #TEST =delta_c_func_my(mesh, 0.5, epsilon, dim) #TODO Just for testing Purpose
    values = delta_cos_func(mesh, 0.5, epsilon, dim)
    #np.testing.assert_array_equal(TEST, values)

    
    import matplotlib.pyplot as plt
    plt.plot(xyz[0], values[nxy//2, nxy//2, :], label="1")
    #plt.plot(xyz[0], TEST[nxy//2, nxy//2, :], linestyle="--", label="h")
    plt.set_xlabel("z ; xy(0.5)")
    plt.set_ylabel("h")
    plt.show()






#%%
if __name__ == '__main__':
    from cdf_project.data_etl.etl_utils import gen_x_nnunif_mesh, gen_nonuniform_segment_delta
    delta_vect, delta_space = gen_nonuniform_segment_delta(125)
    mesh =gen_x_nnunif_mesh(delta_vect)
    x = array_from_mesh(1, mesh)


    import matplotlib.pyplot as plt
    plt.plot(x[0], three_time_phi[1][-1])
    plt.show()


