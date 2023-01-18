# -*- coding: utf-8 -*-
"""
Version 1, Extraction Transformation Loading UTILS
Read setup, generete and write / or read fileds
Genereate and write single (1) solution
@author: enzo
"""
#%% import pakages
import gstools as gs
import fipy as fi

from cdf_project.custom_funcs import  get_logger, benchmark
global logger
logger = get_logger()


#%% 1.  create a mesh
def order_mesh_old(seq, size):
    import numpy as np
    seen = set()
    seen_add = seen.add
    list_out =  [x for x in seq if not (x in seen or seen_add(x))]
    assert len(list_out) ==size
    return np.array(list_out)

def reorder_mesh(seq, size):
    import numpy as np
    array_out = np.unique(seq)
    assert len(array_out) ==size
    return array_out



def array_from_mesh(dim, mesh, lens=None):
    array = []
    if lens == None:
        nx = mesh.nx
        ny = mesh.nx #TODO
        nz = mesh.nx #TODO
    else:
        nx = lens
        ny = lens
        nz = lens

    if dim==1:
        array = [mesh.x.value]
    elif dim == 2:
        xy = mesh.cellCenters.value
        x = reorder_mesh(xy[0],nx)
        y = reorder_mesh(xy[1],ny)
        array = [x,y]
    elif dim ==3:
        xyz = mesh.cellCenters.value
        x = reorder_mesh(xyz[0], nx)
        y = reorder_mesh(xyz[1], ny)
        z = reorder_mesh(xyz[2], nz)
        array = [x,y,z]
    return array






def gen_x_mesh(nx, dx, start = 0.0): #101   # 0.1
    mesh = fi.Grid1D(nx=nx, dx= dx)
    mesh = (mesh + start)
    #x = mesh.cellCenters[0].value  #x = mesh.x.value
    return mesh#, [x]

def gen_x_nnunif_mesh(dxvec, start = 0.0):
    mesh = fi.Grid1D(dx= dxvec)
    mesh = (mesh + start)
    #x = mesh.cellCenters[0].value
    return mesh#, [x]


def gen_xy_mesh(nxy, dxy, start = [[0.0],[0.0]]): #nxy = 101 , dxy =.01
    mesh = fi.Grid2D(nx=nxy, dx= dxy, ny=nxy, dy= dxy)
    mesh = (mesh + start)
    #xy = mesh.cellCenters.value
    #x = order_mesh(xy[0], mesh.nx)
    #y = order_mesh(xy[1], mesh.ny)
    return mesh#, [x, y]

def gen_xy_nnunif_mesh(dxvec, start = [[0.0],[0.0]]):
    mesh = fi.Grid2D(dx= dxvec, dy= dxvec)
    mesh = (mesh + start)
    # xy = mesh.cellCenters.value
    # x = order_mesh(xy[0], len(dxvec))
    # y = order_mesh(xy[1], len(dxvec))
    return mesh#, [x, y]

def gen_xyz_mesh(nxyz , dxyz, start = [[-1.0],[-1.0],[-1.0]]): #nxyz = 101 , dxyz =.01
    logger.debug(f'generating 3D mesh...')
    mesh = fi.Grid3D(dx = dxyz, dy = dxyz, dz = dxyz, \
        nx = nxyz, ny = nxyz, nz = nxyz)
    mesh = (mesh + start)
    logger.debug(f'generated 3D mesh...')
    # xyz = mesh.cellCenters.value
    # x = order_mesh(xyz[0], mesh.nx)
    # y = order_mesh(xyz[1], mesh.ny)
    # z = order_mesh(xyz[2], mesh.nz)
    return mesh#, [x, y, z]

def gen_xyz_nnunif_mesh(dxvec, start = [[-1.0],[-1.0],[-1.0]]):
    mesh = fi.Grid3D(dx = dxvec, dy = dxvec, dz = dxvec)
    mesh = (mesh + start)
    # xyz = mesh.cellCenters.value
    # x = order_mesh(xyz[0], mesh.nx)
    # y = order_mesh(xyz[1], mesh.ny)
    # z = order_mesh(xyz[2], mesh.nz)
    return mesh#, [x, y, z]




#%% 1.a  establish delta

def gen_uniform_delta(number_of_cells, length = 1.0):
    import math
    #digits = int(math.log10(number_of_cells))+1
    delta_x= length/number_of_cells #round(1/number_of_cells, digits)
    logger.info(f'Uniform: number of cells = {number_of_cells}, delta_x={delta_x} ')
    if delta_x < 4*10**-3:
        logger.debug(f'delta_space = {delta_x} < 4*10**-3')
    else:
        logger.critical(f'delta_space = {delta_x} > 4*10**-3')
    #assert(delta_x < 4*10**-3)
    return delta_x

def gen_nonuniform_central_segment_delta(number_of_cells, length = 1.0):
    assert(length > 0)
    L1, L2, L3 = length*0.25, length*0.50, length*0.25
    assert(L1+L2+L3 == length)
    n1, n3 = int(0.10*number_of_cells), int(0.10*number_of_cells)
    n2 = number_of_cells - n1 - n3
    assert(n2>0 and n1>0 and n3>0)
    dx1, dx2, dx3 = L1/n1, L2/n2, L3/n3
    delta_space=min(dx1, dx2, dx3)
    delta_space_max=max(dx1, dx2, dx3)
    logger.info(f'Segment non-uni: number of cells = {number_of_cells}, delta_min={delta_space},delta_max={delta_space_max}')
    if delta_space < 4*10**-3:
        logger.debug(f'delta_space = {delta_space} > 4*10**-3')
    else:
        logger.critical(f'delta_space = {delta_space} > 4*10**-3')
    #assert(delta_space < 4*10**-3)
    delta_vect = [dx1]*n1 + [dx2]*n2 + [dx3]*n3
    return delta_vect, delta_space


def gen_nonuniform_segment_delta(number_of_cells, length = 1.0):
    assert(length > 0)
    L1, L2, L3 = length*0.43, length*0.34, length*0.23
    assert(L1+L2+L3 == length)
    n1, n3 = int(0.07*number_of_cells), int(0.25*number_of_cells)
    n2 = number_of_cells - n1 - n3
    assert(n2>0 and n1>0 and n3>0)
    dx1, dx2, dx3 = L1/n1, L2/n2, L3/n3
    delta_space=min(dx1, dx2, dx3)
    delta_space_max=max(dx1, dx2, dx3)
    logger.info(f'Segment non-uni: number of cells = {number_of_cells}, delta_min={delta_space},delta_max={delta_space_max}')
    if delta_space < 4*10**-3:
        logger.debug(f'delta_space = {delta_space} > 4*10**-3')
    else:
        logger.critical(f'delta_space = {delta_space} > 4*10**-3')
    #assert(delta_space < 4*10**-3)
    delta_vect = [dx1]*n1 + [dx2]*n2 + [dx3]*n3
    return delta_vect, delta_space

def gen_nonuniform_chebyshev_delta(number_of_cells, length = 1.0 ):
    import numpy as np
    assert(length > 0)
    xmin=0.0 # it is infuneltal since there is a diff
    xmax=length
    # This function calculates the n Chebyshev points
    number_of_cells=number_of_cells+1 # one lost becaouse of diff
    ns = np.arange(1,number_of_cells+1)
    x = np.cos((2*ns-1)*np.pi/(2*number_of_cells))
    inverse_verc = (xmin+xmax)/2 + (xmax-xmin)*x/2 #(xmin+xmax)/2 dose note infulece dvec
    dxvec = inverse_verc[:-1]-inverse_verc[1:]
    delta_space=min(dxvec)
    delta_space_max=max(dxvec)
    logger.info(f'Chebyshev: number of cells = {number_of_cells-1}, delta_min={delta_space},delta_max={delta_space_max}')
    assert(delta_space>0.0)
    return dxvec, delta_space


def gen_nonuniform_inverse_chebyshev_delta(number_of_cells, length = 1.0 ):
    dxvec, delta_space = gen_nonuniform_chebyshev_delta(number_of_cells, length)
    #TODO invert dxvec
    #https://www.reddit.com/r/learnpython/comments/gma0is/splitting_an_array_into_two_and_swapping_the/
    
    assert (number_of_cells == len (dxvec))
    #n=number_of_cells//2
    #dxvec[:n],dxvec[n:] = dxvec[n:].copy(), dxvec[:n].copy()  # split_swap_vec
    return dxvec, delta_space




#%% 1.  Field generating with multi processor!

def initializer(input_mesh,ens_no, dims):
    import numpy as np
    sizes = ()
    for array in input_mesh:
        sizes += array.shape
    logger.info(f'sizes={sizes}')
    sizes = (ens_no,*sizes[-dims:])
    fields = np.zeros(sizes)#3D
    logger.info(f'fieldsShape={fields.shape}' )
    return fields

def init_pool(fields):
    global fields_global
    fields_global = fields


def get_result_field(result, start_ensamble):
    global fields_global
    if not result[0]%5: 
        logger.debug(f'merging fileds; working on {result[0]-start_ensamble} i.e. {result[0]}-{start_ensamble}')
    fields_global[result[0]-start_ensamble] = result[1]

def gen_field(i, input_mesh, input_dict, seed, log_normal = True):
    dims = input_dict['dimensions']
    corrx = input_dict['correlation_coeff']
    variance = input_dict['variance']#0.84**2 = 0.7 OR  0.001
    model = gs.Gaussian(dim=dims, var=variance, len_scale= corrx) #VARiance corrx*(lx/ppm)
    srf = gs.SRF(model, mean=0.)
    if not i%50: logger.debug(f'i{i}, seed={seed}')
    if not i%63: logger.info(f'i{i}, seed={seed}')
    srf.structured(input_mesh, mesh_type='structured', seed=seed) #un? structured
    if log_normal:
        gs.transform.normal_to_lognormal(srf)
    return i, srf.field

@benchmark
def multi(input_mesh, input_dict, start_ensamble, ens_no, seeds, log_normal):
    import multiprocessing as mp
    from functools import partial
    global fields_global
    fields_global = initializer(input_mesh, ens_no, input_dict['dimensions'])
    pool = mp.Pool(mp.cpu_count()-1, init_pool(fields_global))
    new_callback_function = partial(get_result_field, start_ensamble=start_ensamble)
    for i in range(start_ensamble, start_ensamble+ens_no):
        pool.apply_async(gen_field, 
        args=(i, input_mesh, input_dict, seeds[i], log_normal), 
        callback=new_callback_function)
    pool.close()
    pool.join()
    return fields_global

@benchmark
def single(input_mesh, input_dict, start_ensamble, ens_no, seeds, log_normal):
    global fields_global
    fields_global = initializer(input_mesh, ens_no, input_dict['dimension'])
    for i in range(start_ensamble, start_ensamble+ens_no):
        get_result_field(gen_field(i, input_mesh, input_dict, seeds[i], log_normal), start_ensamble)
    return fields_global

def ensamble_field(input_mesh, input_dict, ens_no = 1000, start_ensamble=0, log_normal=True): #0.25 0.01   correlation length vector (m)
    #input_mesh = [x, y] is a list of the same lenght as dims
    assert len(input_mesh) == input_dict['dimensions']

    import numpy as np
    from gstools.random import MasterRNG
    seed  = MasterRNG(20170519)
    seeds = np.zeros(start_ensamble+ens_no, np.int64) #2D
    for i in range(0, start_ensamble+ens_no):
        seeds[i]= seed()
    logger.debug(f'generated {seeds.size} size seeds')

    #rng = range(start_ensamble, start_ensamble+ens_no)
    all_fields = multi(input_mesh, input_dict, start_ensamble, ens_no, seeds, log_normal)
    #all_fields = single(input_mesh, input_dict, start_ensamble, ens_no, seeds, log_normal)
    return_dict={}
    return_dict['fields_array'] = all_fields
    return_dict['seeds_array'] = seeds[start_ensamble:(start_ensamble+ens_no)]
    assert len(return_dict['fields_array']) == len(return_dict['seeds_array'])
    return return_dict

# def ensamble_single_field(input_mesh, corrx=0.25, dims= 1, ens_no = 1000): #0.25 0.01   correlation length vector (m)
#     #input_mesh = [x, y] is a list of the same lenght as dims
#     assert len(input_mesh) == dims
#     #%% structured field
#     model = gs.Gaussian(dim=dims, var=0.84**2, len_scale= corrx) #VARiance corrx*(lx/ppm)
#     #model = Stab(dim=1, var=0.84**2, len_scale= corrx*ppm) #VARiance 
#     srf = gs.SRF(model, mean=0.)

#     #%% Ensemble of Fields
#     fields = []

#     from gstools.random import MasterRNG
#     seed = MasterRNG(20170519)
#     for i in range(ens_no):
#         srf.structured(input_mesh, mesh_type='structured', seed=seed()) #un? structured
#         gs.transform.normal_to_lognormal(srf)
#         fields.append(srf.field)

#     return fields



# %% 2. wrapper function

default_dict={
'dimensions' : 1,
'correlation_coeff': 0.01,
'variance':0.7
}

def gen_fields(dims= 1, number_of_cells=128, input_dict=default_dict, ens_num=1000, enstart= 0, log_normal=False,  non_uniform='uniform', interval_len = 1.0, interval_start= 0.0 ):
    assert isinstance(dims, int)
    assert isinstance(number_of_cells, int)
    assert dims <4  and dims>0
    assert (dims == input_dict['dimensions'])

    if non_uniform == 'nonuniform_segment':
        #nonuniform_segment
        delta_vect = gen_nonuniform_segment_delta(number_of_cells, interval_len)
        if dims == 1: 
            mesh =gen_x_nnunif_mesh(delta_vect, interval_start)
        elif dims == 2:
            mesh = gen_xy_nnunif_mesh(delta_vect, interval_start)
        elif dims == 3:
            mesh = gen_xyz_nnunif_mesh(delta_vect, interval_start)
    elif non_uniform == 'uniform':
        delta_space = gen_uniform_delta(number_of_cells, interval_len)
        if dims == 1:
            mesh =gen_x_mesh(number_of_cells, delta_space, interval_start)
        elif dims == 2:
            mesh = gen_xy_mesh(number_of_cells, delta_space, interval_start)
        elif dims == 3:
            mesh = gen_xyz_mesh(number_of_cells, delta_space, interval_start)
    elif non_uniform == 'nonuniform_chebyshev':
        delta_vect, delta_space = gen_nonuniform_chebyshev_delta(number_of_cells, interval_len)
        if dims == 1: 
            mesh =gen_x_nnunif_mesh(delta_vect, interval_start)
        elif dims == 2:
            mesh = gen_xy_nnunif_mesh(delta_vect, interval_start)
        elif dims == 3:
            mesh = gen_xyz_nnunif_mesh(delta_vect, interval_start)
    else:
        raise NameError('Not a valid mesh type')
    
    d_vector_mesh = array_from_mesh(dims, mesh)
    fields = ensamble_field(d_vector_mesh, input_dict, ens_no=ens_num,start_ensamble=enstart, log_normal=log_normal)
    return_dict={}
    return_dict['mesh'] = mesh
    return_dict['fields_array'] = fields['fields_array']
    return_dict['seeds_array'] = fields['seeds_array'] 

    return return_dict





#%% Test single mesh extraction


if __name__ == '__main__':
    import matplotlib.pyplot as pt
    import numpy as np
    num_of_cells = 300
    delta_space = gen_uniform_delta(num_of_cells, length=2)    

    from time import time
    start = time()
    #code here
    mesh =gen_xyz_mesh(num_of_cells, delta_space, start=-1)
    print(f'Time taken to run: {time() - start} seconds')

    xyz = mesh.cellCenters.value
    start = time()
    x1 = reorder_mesh(xyz[0], mesh.nx)
    print(f'Time taken to run order_mesh: {time() - start} seconds')
    
    start = time()
    x2 =np.unique(xyz[0])
    print(f'Time taken to run unique: {time() - start} seconds')

    start = time()
    x3 = list(dict.fromkeys(xyz[0]))
    print(f'Time taken to run list(dict: {time() - start} seconds')

    #np.testing.assert_array_equal(x1, x2)
    np.testing.assert_array_equal(x2, x3)



#%% Test single mesh extraction

if __name__ == '__main__':
    import matplotlib.pyplot as pt
    mesh = fi.Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)

    print(fi.numerix.nonzero(mesh.facesTop)[0])

    xyz = mesh.cellCenters.value
    x = reorder_mesh(xyz[0], mesh.nx)
    y = reorder_mesh(xyz[1], mesh.ny)
    z = reorder_mesh(xyz[2], mesh.nz)


    #np.unique(xyz[0])


#%% Test 0 mesh generations

if __name__ == '__main__':
    import matplotlib.pyplot as pt
    mesh = gen_x_mesh(nx = 505 , dx =.005) #nx = 1001 , dx =.001
    x = array_from_mesh(1, mesh)
    #fields = ensamble_field([x], ens_no = 4, corrx=0.01)
    fields = ensamble_field(x, corrx=0.01, ens_no = 4)['fields_list']
    import numpy as np
    fig, ax = pt.subplots(2, 2, sharex=True, sharey=True)
    ax = ax.flatten()
    for i in range(len(fields)):
        ax[i].plot(fields[i].T)
        #ax[i].set_xticks(np.arange(0, len(x), step=200))
        print(i)
    pt.show()

#%% Test 1 

if __name__ == '__main__':
    number_of_cells = 10    
    delta_vect, delta_space  = gen_nonuniform_chebyshev_delta(number_of_cells, 2.0)
    print(delta_vect)
    mesh=gen_x_nnunif_mesh(delta_vect, start=0.0)
    #print (x)

    delta_vect, delta_space  =  gen_nonuniform_inverse_chebyshev_delta(number_of_cells, 2.0)
    print(delta_vect)
    mesh=gen_x_nnunif_mesh(delta_vect, start=0.0)
    #print (x)



#%% Test 1 
#TODO to debug!!!!
if __name__ == '__main__':
    number_of_cells = 50
    delta_vect, delta_space = gen_nonuniform_segment_delta(number_of_cells)
    mesh=gen_x_nnunif_mesh(delta_vect)
    x = array_from_mesh(1, mesh)
    print (x)

    delta_vect, delta_space = gen_nonuniform_segment_delta(number_of_cells, length=2)
    mesh=gen_x_nnunif_mesh(delta_vect, start=-1)
    x = array_from_mesh(1, mesh)
    print (x)

#%% Test 1 

if __name__ == '__main__':
    number_of_cells = 10
    delta_space = gen_uniform_delta(number_of_cells)
    mesh =gen_x_mesh(number_of_cells, delta_space)
    x = array_from_mesh(1, mesh)
    print (x)
    delta_space = gen_uniform_delta(number_of_cells, length=2)
    mesh =gen_x_mesh(number_of_cells, delta_space, start=-1)
    x = array_from_mesh(1, mesh)
    print (x)

    delta_space = gen_uniform_delta(number_of_cells, length=2)
    mesh =gen_x_mesh(number_of_cells, delta_space, start=-1)
    x = array_from_mesh(1, mesh)
    print (x)

    mesh  = gen_xy_mesh(number_of_cells, delta_space, start=-2)
    d_vector_mesh = array_from_mesh(2, mesh)
    print (d_vector_mesh)

    mesh = gen_xy_mesh(number_of_cells, delta_space, start=[[1], [4]])
    d_vector_mesh = array_from_mesh(2, mesh)
    print (d_vector_mesh)


# %% Test 2D

if __name__ == '__main__':
    import matplotlib.pyplot as pt
    mesh = gen_xy_mesh(nxy = 505 , dxy =.005) #nx = 1001 , dx =.001
    xy = array_from_mesh(2, mesh)
    fields = ensamble_field(xy, corrx=0.01, dims=2, ens_no = 4)['fields_list']

    fig, ax = pt.subplots(2, 2, sharex=True, sharey=True)
    ax = ax.flatten()
    for i in range(len(fields)):
        print(i)
        #print(field[i])
        ax[i].imshow(fields[i].T, origin="lower")
    pt.show()


# %% Test 2D in 3D
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    mesh = gen_xy_mesh(nxy = 505 , dxy =.005) #nx = 1001 , dx =.001
    xy = array_from_mesh(2, mesh)
    fields = ensamble_field(xy, corrx=0.01, dims=2, ens_no = 2)['fields_list']

    # Twice as wide as it is tall.
    fig = plt.figure(figsize=plt.figaspect(0.5))
    X, Y = np.meshgrid(xy[0], xy[1])

    #---- First subplot
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X, Y, fields[0].T, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=10)

    #---- Second subplot
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_wireframe(X, Y, fields[0], rstride=10, cstride=10)
    plt.show()




# %%
