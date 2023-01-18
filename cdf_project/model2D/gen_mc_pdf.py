#%% Import
import numpy as np
import matplotlib.pyplot as pt

from utils_field import gen_x_mesh, ensamble_field
from utils_pde import pde_f_time_fun, verify_param
from utils_pdf import C_p


def get_result_time(result):
    global allMC_phi_tt
    global dt_array
    allMC_phi_tt[result[0]] = result[1]
    dt_array[result[0]] = result[2]
 


#%% Mesh and field
dx =.002  #001  #nx = 101 , dx =.01
nx = 501 # 501 .002
mesh, x = gen_x_mesh(nx = nx , dx =dx) #nx = 1001 , dx =.001

MC_runs= 100 #1000
correl=0.001 # 0.0001  25  001
correl_text='Corr' #  NonCorr  Corr
cr_txt='c_0filed' #  nc c



fields = ensamble_field(x, ens_no = MC_runs, corrx=correl)['fields_list']


# %%  pre assig parameters + convergence condition by Courant

print(f'x = {len(x)}')
print(f'field = {len(fields[0])}')
assert len(x) == len(fields[0])
size_mesh_x = len(x)
mc_runs = len(fields)
assert mc_runs == MC_runs

D = 1./100.  #1./100. 1.      0.75 #0.1


#TODO Estimate overall Courant number, mean value and fix min steps as mean number!

min_time_steps = 150



# %% Final time required

totalElapsedTime=0.1

allMC_phi_tt = np.zeros((len(fields), 3+1, size_mesh_x))
dt_array = np.zeros((len(fields), 1))


import multiprocessing as mp
import time
ts = time.time()
pool = mp.Pool(mp.cpu_count()-1)
for i in range(mc_runs):
    pool.apply_async(pde_f_time_fun, 
    args=(i, D, mesh, fields[i],totalElapsedTime,  False), 
    callback=get_result_time)
pool.close()
pool.join()
print('Time in parallel:', time.time() - ts)



#%% for Seq
totalElapsedTime=0.1
allMC_phi_tt = np.zeros((len(fields), 3+1, size_mesh_x))
dt_array = np.zeros((len(fields), 1))

import time
ts = time.time()
for i in range(0, mc_runs):
    get_result_time(pde_f_time_fun(i, D, mesh, fields[i], totalElapsedTime,\
        verbose=True))  #False  True
print('Time in serial:', time.time() - ts)


     
#%% 7. control

print(allMC_phi_tt.shape)

_ , Time, Space = allMC_phi_tt.shape

m= Time-1    #150  # 0- 150
n= Space//2   #  50  # 0-101


#print(sample)

#%% 7. Save

import pickle
file_name=f'allMC_phi_tt_{mc_runs}_{cr_txt}_{nx}.plk'

#%% 7. Save
#with open(file_name, 'wb') as f:
#    pickle.dump(allMC_phi_tt, f)

#%% 7. Open
#file_name=f'allMC_phi_tt_{mc_runs}_nc_{nx}.plk'
#file_name=f'allMC_phi_tt_{mc_runs}_cr_{nx}.plk'

with open(file_name, 'rb') as f:
    allMC_phi_tt = pickle.load(f)


#%% 7. search max n index!
m = 0
sample = allMC_phi_tt[:, m, :]
print(sample.shape)
max_positions = np.argmax(sample, axis=1)
print(max_positions.shape)
print(f'time = {m}  mean={np.mean(max_positions)}, std={np.std(max_positions)}')

m = 1
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

m = 3
sample = allMC_phi_tt[:, m, :]
print(sample.shape)
max_positions = np.argmax(sample, axis=1)
print(max_positions.shape)
print(f'time = {m}  mean={np.mean(max_positions)}, std={np.std(max_positions)}')




def m_out(m):
    if m==2: return 'Final'
    elif m==1: return 'mid'
    else : return 'start'


#%% 7. CDF plot

import matplotlib.pyplot as plt

n_0=250
n_1=250
n_2=250
n_3=250



# plot the sorted data:
fig = plt.figure(figsize=(10, 5))
fig.suptitle(f'D={D}, totalTime={totalElapsedTime}, dx={dx}, {correl_text}', fontsize=9)




ax1 = fig.add_subplot(121)
data_sorted, p = C_p(allMC_phi_tt, 1, n_1) #n
ax1.plot(data_sorted, p, label=f't≈{totalElapsedTime*0.1:.2f};x={n_1}')
data_sorted, p = C_p(allMC_phi_tt, 2, n_2)
ax1.plot(data_sorted, p, label=f't≈{totalElapsedTime*0.5:.2f};x={n_2}')
data_sorted, p = C_p(allMC_phi_tt, 3, n_3)
ax1.plot(data_sorted, p, label=f't≈{totalElapsedTime:.2f};x={n_3}')

ax1.set_xlabel('$C$') #SCALA 0- 1 E NO 0.2 - 0.24, PRENDERE IL C MAX OVVERO A 0.5
ax1.set_ylabel('$p$')
ax1.set_title(f'relative')
ax1.set_ylim(0,1)
ax1.legend()





ax2 = fig.add_subplot(122)
data_sorted, p = C_p(allMC_phi_tt, 1, n_1) #n
Cmin = data_sorted.min()
Cmax = data_sorted.max()
data_sorted=data_sorted/Cmax
ax2.plot(data_sorted, p, label=f'Cmin={Cmin:.0f} Cmax={Cmax:.2f} t≈{totalElapsedTime*0.1:.2f};x={n_1*dx:.2f}')

data_sorted, p = C_p(allMC_phi_tt, 2, n_2)
Cmin = data_sorted.min()
Cmax = data_sorted.max()
data_sorted=data_sorted/Cmax
ax2.plot(data_sorted, p, label=f'Cmin={Cmin:.0f} Cmax={Cmax:.2f} t≈{totalElapsedTime*0.5:.2f};x={n_2*dx:.2f}')


data_sorted, p = C_p(allMC_phi_tt, 3, n_3)
Cmin = data_sorted.min()
Cmax = data_sorted.max()
data_sorted=data_sorted/Cmax
ax2.plot(data_sorted, p, label=f'Cmin={Cmin:.0f} Cmax={Cmax:.2f} t≈{totalElapsedTime:.2f};x={n_3*dx:.2f}')


ax2.set_xlabel('$C$') #SCALA 0- 1 E NO 0.2 - 0.24, PRENDERE IL C MAX OVVERO A 0.5
ax2.set_ylabel('$p$')
ax2.set_title(f'normalized')
ax2.set_ylim(0,1)
ax2.set_xticks(np.arange(0, 1.1, step=1/10))
ax2.legend()


plt.savefig(f'CDF_{cr_txt}_{n_1}_{n_3}.png', dpi=300)
plt.show()




#%% 7. cumulative:  1.mean - 2.std - 3.   plot
#https://stackoverflow.com/questions/53005146/calculating-the-cumulative-mean-in-python

probe = allMC_phi_tt[1:, 1, nx//2] # almeno 2


fig, axs = plt.subplots(2, 2, figsize=(10, 10))

eq_formula = fr'$\alpha > \beta {D}$'

fig.suptitle(f'{eq_formula}, dx= {dx}, fulltime≈ {totalElapsedTime:.2f} \n {correl_text}  - over MC runs:{mc_runs}')


#mean
cum_summ_dr = np.cumsum(probe)
N = np.arange(2,mc_runs+1)
N1 = np.arange(1,mc_runs)
mean_csmdr = cum_summ_dr / N
axs[0, 0].plot(N, mean_csmdr, 'tab:gray')
axs[0, 0].set_title('Mean')

#variance - second order
var=probe- mean_csmdr
var_res =  np.cumsum(var**2)/ N1
std_res = np.sqrt(var_res)
axs[0, 1].plot(x, var_res, 'tab:red')
axs[0, 1].set_title('Variance')

#skewness - third order
skew_res =  np.cumsum((var/std_res)**3)/ N
axs[1, 0].plot(x, skew_res, 'tab:green')
axs[1, 0].set_title('Skewness')


#kurtosis - 4th-order
kurt_res =  np.cumsum((var/std_res)**4)/ N - 3
axs[1, 1].plot(x, kurt_res, 'tab:orange')
axs[1, 1].set_title('Kurtosis')

plt.tight_layout()
plt.savefig(f'moments_{cr_txt}.jpg', dpi=300)



# %% Plot function

def m_size(m, totalElapsedTime):
    if m==2: return totalElapsedTime*0.5
    elif m==1: return totalElapsedTime*0.1
    elif m ==0: return 0
    elif m ==3: return totalElapsedTime
    else : raise AssertionError()


fig, axs = plt.subplots(2, 2, sharex=False) #, sharex=True
fig.suptitle(f'D={D}, dx= {dx}, fulltime≈ {totalElapsedTime:.2f} - {correl_text}  - over MC runs:{mc_runs}')
plt.tight_layout()


m=0
for ax in axs.flat:

    for i in range(mc_runs):
        y=allMC_phi_tt[i, m, :]
        ax.plot(np.linspace(0,1,nx), y, marker = '+', linestyle = '', alpha=0.1 )

    ax.set_xlabel(f'time ≈ {m_size(m, totalElapsedTime):.2f}',fontsize='small' )
    ax.tick_params(axis='both', which='major', labelsize=7)
    #axs[m].set_ylabel(f'y')
    ax.set_xlim([0,1])
    print(m)
    m+=1



print('render....')
plt.savefig(f'functs_{cr_txt}_times.jpg', dpi=300)
print('end')




# %%
