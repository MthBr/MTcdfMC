#%% Import
import gstools as gs
import fipy as fi
import numpy as np

verbose = True #False

#%% def internal funcrtion

def values_4_delta(vector = [0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085,
       0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175,
       0.185, 0.195, 0.205, 0.215, 0.225, 0.235, 0.245, 0.255, 0.265,
       0.275, 0.285, 0.295, 0.305, 0.315, 0.325, 0.335, 0.345, 0.355,
       0.365, 0.375, 0.385, 0.395, 0.405, 0.415, 0.425, 0.435, 0.445,
       0.455, 0.465, 0.475, 0.485, 0.495, 0.505, 0.515, 0.525, 0.535,
       0.545, 0.555, 0.565, 0.575, 0.585, 0.595, 0.605, 0.615, 0.625,
       0.635, 0.645, 0.655, 0.665, 0.675, 0.685, 0.695, 0.705, 0.715,
       0.725, 0.735, 0.745, 0.755, 0.765, 0.775, 0.785, 0.795, 0.805,
       0.815, 0.825, 0.835, 0.845, 0.855, 0.865, 0.875, 0.885, 0.895,
       0.905, 0.915, 0.925, 0.935, 0.945, 0.955, 0.965, 0.975, 0.985,
       0.995, 1.005]):
    
    middleIndex = (len(vector) - 1)//2
    middle_value =  vector[middleIndex]
    #if verbose : 
    #print(f'midd= {middleIndex}, mid value = {middle_value} \n')
    delta = vector[middleIndex+1]- middle_value
    assert delta>0
    xa = vector[middleIndex] - delta*0.7
    xb = vector[middleIndex] + delta*0.7
    assert xa < middle_value
    assert middle_value < xb
    return (xa, xb, delta)


#print(values_4_delta(vector = [ 0.485, 0.495, 0.505, 0.515, 0.525, 0.535]))

#print(values_4_delta(vector = [0.495, 0.505, 0.515, 0.525, 0.535]))


def get_result_test(result):
    global allMC_phi_tt_local
    allMC_phi_tt_local.append(result[1])


def get_result(result):
    global allMC_phi_tt_local
    allMC_phi_tt_local[result[0]] = result[1]



#%% def functions
def pde_moive_function(i, D, mesh, field, timeStepDuration, steps=150, verbose=False):
    #%% 4.  set Param
    x = mesh.cellCenters[0].value #TODO togliere x??? SERV SOLO per phi inzi

    U= fi.CellVariable(mesh=mesh, rank=1)
    U[:] =  field #- (1-x) 

    #%% 2.  Initial conditions
    #phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=3.)
    ###u = fi.CellVariable(name="u",  hasOld=True, mesh=mesh, value=0.)
    ###phi.setValue(1., where=(0.25 < x) & (.75>x))

    #Initial Impulse condition
    phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    #phi.setValue(1., where=(.499<x) & (.509>x))
    xa, xb, dx = values_4_delta(x) #TODO spostare calcolo fuori da pdeFunc, perche stesso per onig
    # un subset di paramet che p func=phi.. (kwarg* xa, xb)
    phi.setValue(1., where=(xa<x) & (xb>x)) #??? per la x???



    #%% 3. boundary conditions
    #phi.constrain(0., mesh.facesLeft)
    #phi.constrain(0., mesh.facesRight)

    BCs = (fi.FixedFlux(faces=mesh.facesRight, value=0.),  
            fi.FixedFlux(faces=mesh.facesLeft, value=0.))
    #fixed flux (Neumann condition)

    #%% Equation
    equ =  (fi.TransientTerm(var=phi) ==
            fi.DiffusionTerm(coeff=D, var=phi)
            - fi.ConvectionTerm(coeff=U, var=phi))


    #%% 7. Time setup & solver
    #timeStepDuration = 1.e-4#0.0001
    #steps = 150 # 2000
    
    mySolver = fi.LinearPCGSolver(iterations=1234, tolerance=5e-6)
    #mySolver = fi.LinearGMRESSolver(iterations=1000,tolerance=1.e-6)

    #https://fipy.nist.narkive.com/wIKBXdLD/residuals-in-fipy

    desiredResidual = 1e-9
    elapsedTime = 0
    tryMax=10

    if verbose:
         vi = fi.Viewer((phi), limits={'xmin': 0.0, 'xmax': 1.0})
         viewerMPL = fi.Matplotlib1DViewer(phi,\
             limits={'xmin': 0.0, 'xmax': 1.0},\
             title=f"Delta")

    #%% 7. Solving PDE
    phi_tt = np.zeros((steps+1, len(phi.value)))
    phi_tt[0]  = phi.value.copy()
    #while elapsedTime < totalElapsedTime:
    for step in range(steps):
        phi.updateOld()

        residual=10
        try_count=0
        while residual> desiredResidual and try_count < tryMax:
            residual = equ.sweep(var=phi, dt=timeStepDuration, boundaryConditions = BCs)
            try_count +=1
            if verbose:  print(f'residual = {residual},try_count = {try_count} ')
        assert(try_count < tryMax-1)
        
        elapsedTime += timeStepDuration
        #print(f'elapsedTime = {elapsedTime}')
        if verbose: 
            vi.plot()
            viewerMPL.plot()
            print(f'time_step = {step}')
        phi_tt[step+1] = phi.value.copy()
        #phi_tt.append(phi.value.copy())
    
    if not i%20:
        print(i)
    return (i, phi_tt)





def pde_f_time_fun(i, D, mesh, field, totalElapsedTime=0.1, verbose=False):
    #%% 4.  set Param
    x = mesh.cellCenters[0].value #TODO togliere x??? SERV SOLO per phi inzi

    U= fi.CellVariable(mesh=mesh, rank=1)
    U[:] =  field #- (1-x) 

    #F= fi.Variable(mesh=mesh, rank=1)
    #F[:] = 1* field #- (1-x) 


    #%% 2.  Initial conditions
    #phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=3.)
    ###u = fi.CellVariable(name="u",  hasOld=True, mesh=mesh, value=0.)
    ###phi.setValue(1., where=(0.25 < x) & (.75>x))

    #Initial Impulse condition
    phi = fi.CellVariable(name="phi", mesh=mesh, hasOld=True,  value=0.)
    #phi.setValue(1., where=(.499<x) & (.509>x))
    xa, xb, dx = values_4_delta(x) #TODO spostare calcolo fuori da pdeFunc, perche stesso per onig
    # un subset di paramet che p func=phi.. (kwarg* xa, xb)
    phi.setValue(1., where=(xa<x) & (xb>x)) #??? per la x???



    #%% 3. boundary conditions
    #phi.constrain(0., mesh.facesLeft)
    #phi.constrain(0., mesh.facesRight)

    BCs = (fi.FixedFlux(faces=mesh.facesRight, value=0.),  
            fi.FixedFlux(faces=mesh.facesLeft, value=0.))
    #fixed flux (Neumann condition)

    #%% Equation
    equ =  (fi.TransientTerm(var=phi) ==
            fi.DiffusionTerm(coeff=D, var=phi) 
            + U
            - fi.ConvectionTerm(coeff=U, var=phi))


    #%% 7. Time setup & solver
    if verbose: print(f'dx={dx}')
    u_max = max(field)
    

    safetyFactor = 0.9
    timestepDuration_old = 0.9*dx/u_max
    timestepDuration = safetyFactor * dx**2 / (2 * u_max)
    print(timestepDuration_old)
    print(timestepDuration)

    time_steps = int(totalElapsedTime//timestepDuration_old)
    #time_steps = max(300, time_steps)
    timeStepDuration = totalElapsedTime/time_steps

    print(f'time_steps = {time_steps}')
    print(f'timeStepDuration = {timeStepDuration}')


    
    mySolver = fi.LinearPCGSolver(iterations=1234, tolerance=5e-6)
    #mySolver = fi.LinearGMRESSolver(iterations=1000,tolerance=1.e-6)

    #https://fipy.nist.narkive.com/wIKBXdLD/residuals-in-fipy

    desiredResidual = 1e-9
    elapsedTime = 0
    tryMax=10

    #%% 7. Solving PDE
    phi_tt = np.zeros((3+1, len(phi.value)))
    phi_tt[0]  = phi.value.copy()


    #while elapsedTime < totalElapsedTime:
    for step in range(time_steps//10):
        phi.updateOld()

        residual=10
        try_count=0
        while residual> desiredResidual and try_count < tryMax:
            residual = equ.sweep(var=phi, dt=timeStepDuration, boundaryConditions = BCs)
            try_count +=1
            #if verbose:  print(f'residual = {residual},try_count = {try_count} ')
        assert(try_count < tryMax-1)
        elapsedTime += timeStepDuration
    phi_tt[1] = phi.value.copy()
    if verbose: print(f'MIN elapsedTime = {elapsedTime}, steps{time_steps//10}')


    #while elapsedTime < totalElapsedTime:
    for step in range(time_steps//10, time_steps//2):
        phi.updateOld()

        residual=10
        try_count=0
        while residual> desiredResidual and try_count < tryMax:
            residual = equ.sweep(var=phi, dt=timeStepDuration, boundaryConditions = BCs)
            try_count +=1
            #if verbose:  print(f'residual = {residual},try_count = {try_count} ')
        assert(try_count < tryMax-1)
        elapsedTime += timeStepDuration
    phi_tt[2] = phi.value.copy()
    if verbose: print(f'MID elapsedTime = {elapsedTime}')
    
    for step in range(time_steps//2):
        phi.updateOld()

        residual=10
        try_count=0
        while residual> desiredResidual and try_count < tryMax:
            residual = equ.sweep(var=phi, dt=timeStepDuration, boundaryConditions = BCs)
            try_count +=1
            #if verbose:  print(f'residual = {residual},try_count = {try_count} ')
        assert(try_count < tryMax-1)
        
        elapsedTime += timeStepDuration
    phi_tt[3] = phi.value.copy()
    if verbose: print(f'FIN elapsedTime = {elapsedTime}')

    if not i%20:
        print(f'{i}, fin_time={elapsedTime}, CouranrNum={u_max*timeStepDuration/dx:.5f}')
        print(f'{i}, minstep={time_steps//10}, midStep={time_steps//2}, ENDStep={time_steps}')
        print(f'{i}, dtnew{timeStepDuration}  dt old {timeStepDuration_old}')
    return (i, phi_tt,timeStepDuration)



def verify_param(fields, timeStepDuration, dx):
    u_max = max(fields, key=tuple).max()
    u_min = min(fields, key=tuple).min()
    u_mean = np.mean(fields[:])
    #Cour = D*
    print(f'CourMAX={u_max*timeStepDuration/dx}, \
        CourMean={u_mean*timeStepDuration/dx}, \
            CourMIN={u_min*timeStepDuration/dx}')




#%% Test


if __name__ == '__main__':
    import numpy as np
    from utils_field import gen_x_mesh, ensamble_field
    dx =.01
    mesh, x = gen_x_mesh(nx = 101 , dx =dx) #nx = 1001 , dx =.001
    fields = ensamble_field(x, ens_no = 4, corrx=0.01)['fields_list']
    assert len(x) == len(fields[0])
    size_mesh_x = len(x)

    time_steps = 150 # 2000
    allMC_phi_tt_local = np.zeros((len(fields), time_steps+1, size_mesh_x))

    timeStepDuration = 55.e-5# 5.e-5 0.0001
    D = 1./100.  #1./100. 1.      0.75 #0.1

    verify_param(fields, timeStepDuration, dx)


    import time
    ts = time.time()
    for i in range(0, len(fields)):
        get_result(pde_function(i, D, mesh, fields[i], timeStepDuration, time_steps, verbose=False))
    print('Time in serial:', time.time() - ts)
    print(allMC_phi_tt_local)



    import multiprocessing as mp
    allMC_phi_tt_local = np.zeros((len(fields), time_steps+1, size_mesh_x))
    ts = time.time()
    pool = mp.Pool(mp.cpu_count()-1)
    for i in range(0, len(fields)):
        pool.apply_async(pde_function, 
        args=(i, D, mesh, fields[i], timeStepDuration, time_steps, False), 
        callback=get_result)
    pool.close()
    pool.join()
    print('Time in parallel:', time.time() - ts)
    print(allMC_phi_tt_local)




# %%





# %%
