#%% Setup solver
import os
os.environ["FIPY_SOLVERS"] = "scipy" #pyamgx petsc scipy no-pysparse trilinos  pysparse pyamg

os.environ["FIPY_VERBOSE_SOLVER"] = "1" # 1:True # Only for TESTING
#print('Verbose ' + os.environ.get('FIPY_VERBOSE_SOLVER'))
os.environ["OMP_NUM_THREADS"]= "1"


#%% Setup parameters

_test_set_cofig = 100
_field_config = 30 #30 10 40
_mesh_config= 5
_dim_set = 2000
_delta = 1


#Defaults! Nones
id= f'Test_00' #TODO unused#
dxyz = None #TODO remove
mc_runs = 1 #TODO delte
delta_type = 'pg5-s25' # l  cos Cubic LL  3*-f_arcs  # 2-l  1-l  1*-l 2-cos 2-Cubic 2-LL
delta_product = 'tensor' # 'radial' 'scaled_tensor' tensor





solver_type = 'AMGx'

if _mesh_config==5:
    grid_type='uniform' #TODO add
    nxyz = 512 #256 #126  #512 1024
    center = 0
    leng = 2
    start = -1
if _mesh_config==1:
    dxyz = None #TODO remove
    nxyz = 1024 #256 #126  #512 1024
    center = 0
    leng = 2
    start = -1
if _mesh_config==2:
    dxyz = None #TODO remove
    nxyz = 256 #256 #126  #512 1024
    center = 0
    leng = 2
    start = -1



if _field_config==30:
    correl = 0.2   # if nCOR= 0.001 = 10**-3
    var=0.01
    correl_sm_txt = 'cor_detrm'
    correl_long_txt = 'Corr-Determ'
    is_log_normal = False

if _field_config==10:
    correl = 0.2  # if nCOR= 0.001 = 10**-3
    var=2
    correl_sm_txt = 'cor'
    correl_long_txt = 'Corr'
    is_log_normal = False

if _field_config==20:
    correl = 0.1  # if nCOR= 0.001 = 10**-3
    var=0.01
    correl_sm_txt = 'noncor_detrm'
    correl_long_txt = 'NonCorr-Determ'
    is_log_normal = False

if _field_config==40:
    correl = 0.1   # if nCOR= 0.001 = 10**-3
    var=3.79
    correl_sm_txt = 'noncorr'
    correl_long_txt = 'NonCorr'
    is_log_normal = False



if _dim_set == 2000: #TODO 2000
    dims = 2

if _dim_set == 3000:
    dims = 3


if _test_set_cofig==200:
    #dims = 2
    static_prob = True
    BCs = 'None' # Dirichlet Neumann None

    D = 0.0
    convect = False
    source = 'delta' #delta   delta Fdelta  eFdelta 
    phi_0 = 'None'

    eq_formula = fr'$ \nabla \cdot (e^{{Y}} \nabla h) = delta^{{{delta_type}}}(x-0), x\in[{start},{start+leng}]^{dims}\subset R^{dims}$'
    eq_formula2 = rf'{BCs}=0;  filed Y: IntegralScale:{correl}, var:{var}; mesh:{nxyz}'



if _test_set_cofig==100:
    #dims = 2
    static_prob = True
    BCs = 'Dirichlet' # Dirichlet Neumann None

    D = 0.0
    convect = False
    source = 'delta' #delta   delta Fdelta  eFdelta 
    phi_0 = 'None'

    eq_formula = fr'$ \nabla \cdot (e^{{Y}} \nabla h) = delta^{{{delta_type}}}(x-0), x\in[{start},{start+leng}]^{dims}\subset R^{dims}$'
    eq_formula2 = rf'{BCs}=0;  filed Y: IntegralScale:{correl}, var:{var}; mesh:{nxyz}'


if _test_set_cofig==700:
    #dims = 2
    static_prob = True
    BCs = 'Neumann' # Dirichlet Neumann None

    D = 1.0
    convect = True
    source = 'eFdelta' #delta   delta Fdelta  eFdelta 
    phi_0 = 'None'
    
    eq_formula = fr'$ \nabla \cdot ({D} \nabla h) + \nabla Y \cdot \nabla h = e^{{-Y}} * delta^{{{delta_type}}}(x-0), x\in[{start},{start+leng}]^{dims}\subset R^{dims}$'
    eq_formula2 = rf'{BCs}=0;  filed Y: IntegralScale:{correl}, var:{var}; mesh:{nxyz}'


if _test_set_cofig==500:
    #dims = 2
    static_prob = True
    BCs = 'Dirichlet' # Dirichlet Neumann None

    D = 1.0
    convect = True
    source = 'eFdelta' #delta   delta Fdelta  eFdelta 
    phi_0 = 'None'
    
    eq_formula = fr'$ \nabla \cdot ({D} \nabla h) + \nabla Y \cdot \nabla h = e^{{-Y}} * delta^{{{delta_type}}}(x-0), x\in[{start},{start+leng}]^{dims}\subset R^{dims}$'
    eq_formula2 = rf'{BCs}=0;  filed Y: IntegralScale:{correl}, var:{var}; mesh:{nxyz}'



id=  f'Test_{_test_set_cofig+_field_config+_mesh_config+_dim_set+_delta}'

# controls
import math
assert isinstance(correl, float)
assert isinstance(mc_runs, int)
assert isinstance(dims, int)
assert isinstance(D, float)
assert isinstance(convect, bool)
assert isinstance(is_log_normal, bool)
assert isinstance(eq_formula, str)
#assert isinstance(text_long, str) #eq_formula2
#assert isinstance(compact_text, str)


#TODO revisit dxyz
#if dxyz!= None: 
#    assert math.isclose(dxyz*nxyz,1.0, rel_tol=10**-3)

# assign
setup_dict = {}

#setup_dict['grid_spacing'] = dxyz
#setup_dict['center_spacing'] = dxyz #used for delta Dirac
setup_dict['id'] = id
setup_dict['number_of_cells'] = nxyz
setup_dict['length'] = leng
setup_dict['start'] = start


setup_dict['correlation_coeff'] = correl
setup_dict['variance'] = var
setup_dict['correl_text'] = correl_sm_txt
setup_dict['log_normal'] = is_log_normal



setup_dict['dimensions'] = dims

setup_dict['D'] = D
setup_dict['convection'] = convect
setup_dict['initial'] = phi_0


setup_dict['delta_type'] = delta_type
setup_dict['delta_product'] = delta_product



setup_dict['source'] = source 
setup_dict['boundaryConditions'] = BCs 

setup_dict['eq_formula'] = eq_formula
setup_dict['eq_formula2'] = eq_formula2  ###eq_formula2    text_long


setup_dict['compactxt_field'] =  f'cor{int(correl*1000)} variance {int(var*1000)}'  #
setup_dict['compactxt_mesh'] =  f'mesh{nxyz}'  #
setup_dict['compactxt_eq'] =  f'D{int(D*10)}_dlt:{delta_type} GMRE{os.environ["FIPY_SOLVERS"]}'  #

setup_dict['solver_type'] =  f'solver_type'  #

setup_dict['filed_hash_name'] = f"dim{dims}_isLN{is_log_normal}_istr{start}_ilen{leng}_cor{int(correl*100)}_var{int(var*100)}_n{nxyz}"
setup_dict['sol_hash_name'] = f"{id}_{dims}d_{setup_dict['compactxt_eq']}_{setup_dict['compactxt_mesh']}"


setup_dict['hash_log_name'] = f'{id}_dim{dims}'







# %%
