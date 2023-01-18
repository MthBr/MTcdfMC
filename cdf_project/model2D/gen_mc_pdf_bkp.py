#%% Import
import numpy as np
import matplotlib.pyplot as pt

from utils_field import gen_x_mesh, ensamble_field
from utils_pde import pde_moive_function, pde_f_time_fun, verify_param
from utils_pdf import C_p


def get_result(result):
    global allMC_phi_tt
    allMC_phi_tt[result[0]] = result[1]


#%% Mesh and field
dx =.002  #001  #nx = 101 , dx =.01
nx = 501 # 501 .002
mesh, x = gen_x_mesh(nx = nx , dx =dx) #nx = 1001 , dx =.001

MC_runs=1000
correl=0.001
fields = ensamble_field(x, ens_no = MC_runs, corrx=correl)['fields_list']


# %%  pre assig parameters + convergence condition by Courant


print(f'x = {len(x)}')
print(f'field = {len(fields[0])}')
assert len(x) == len(fields[0])
size_mesh_x = len(x)
mc_runs = len(fields)
assert mc_runs == MC_runs


# %% Time steps props
time_steps = 1000  #300 # 2000
allMC_phi_tt = np.zeros((len(fields), time_steps+1, size_mesh_x))

timeStepDuration = 13.e-5 #  13.e-5    7.e-5   3.e-5   5.e-5   0.0001
D = 1./100.  #1./100. 1.      0.75 #0.1

verify_param(fields, timeStepDuration, dx)


#%% for Parallel
import multiprocessing as mp
import time
ts = time.time()
pool = mp.Pool(mp.cpu_count()-1)
for i in range(mc_runs):
    pool.apply_async(pde_moive_function, 
    args=(i, D, mesh, fields[i],timeStepDuration, time_steps,  False), 
    callback=get_result)
pool.close()
pool.join()
print('Time in parallel:', time.time() - ts)


#%% for Seq
import time
ts = time.time()
for i in range(0, mc_runs):
    get_result(pde_moive_function(i, D, mesh, fields[i], timeStepDuration, time_steps,\
        verbose=False))  #False  True
print('Time in serial:', time.time() - ts)


     
#%% 7. control

print(allMC_phi_tt.shape)

_ , Time, Space = allMC_phi_tt.shape

m= Time-1    #150  # 0- 150
n= Space//2   #  50  # 0-101


#print(sample)


#%% 7. search max n index!
m = Time-1
sample = allMC_phi_tt[:, m, :]
print(sample.shape)

max_positions = np.argmax(sample, axis=1)

print(max_positions.shape)
print(f'time = {m}  mean={np.mean(max_positions)}, std={np.std(max_positions)}')

m = Time//2
sample = allMC_phi_tt[:, m, :]
print(sample.shape)

max_positions = np.argmax(sample, axis=1)

print(max_positions.shape)
print(f'time = {m}  mean={np.mean(max_positions)}, std={np.std(max_positions)}')


m = 2
sample = allMC_phi_tt[:, m, :]
print(sample.shape)

max_positions = np.argmax(sample, axis=1)

print(max_positions.shape)
print(f'time = {m}  mean={np.mean(max_positions)}, std={np.std(max_positions)}')




#%% 7. CDF plot

import matplotlib.pyplot as plt


# plot the sorted data:
fig = plt.figure()


data_sorted, p = C_p(allMC_phi_tt, Time-1, 541) #n
plt.plot(data_sorted, p, label=f't={Time-1}')
data_sorted, p = C_p(allMC_phi_tt, Time//2, 525)
plt.plot(data_sorted, p, label=f't={Time//2}')
data_sorted, p = C_p(allMC_phi_tt, 2, n)
plt.plot(data_sorted, p, label=f't={2}')

plt.xlabel('$C$') #SCALA 0- 1 E NO 0.2 - 0.24, PRENDERE IL C MAX OVVERO A 0.5
plt.ylabel('$p$')
plt.title(f'Position={n}')

Cmin = allMC_phi_tt[:, :, n].min()
Cmax = allMC_phi_tt[:, :, n].max()
print(f'Cmin={Cmin}, Cmax={Cmax}')
left, right = plt.xlim()  # return the current xlim
plt.xlim((min(left, Cmin), max(right, Cmax)))   # set the xlim to left, right
plt.ylim(0,1)

plt.xticks(np.arange(Cmin, Cmax, step=(Cmax-Cmin)/5))
plt.legend()

#%% 7. cumulative:  1.mean - 2.std - 3.   plot
#https://stackoverflow.com/questions/53005146/calculating-the-cumulative-mean-in-python





# %%
m=0
m=Time//2
m=Time-1



fig= plt.figure()
for i in range(mc_runs):
    y=allMC_phi_tt[i, m, :]
    plt.plot(np.linspace(0,1,nx), y, marker = '+', linestyle = '', alpha=0.1 )

plt.xlabel(f'x; dx={dx}')
plt.ylabel('y over 1000')
plt.title(f'D={D}, time step= {m}, time≈ {m*timeStepDuration:.2f}, dt={timeStepDuration:.0e} NC')
plt.show()


















# %%

import matplotlib.pyplot as plt
plt.plot(np.linspace(0,1,nx), phi.value)
plt.xlim(0,1)
plt.ylim(0,1)
plt.show()













# %%
