#%% Import
import numpy as np
import matplotlib.pyplot as pt

from utils_field import gen_xy_mesh, ensamble_field
from utils_pde import pde_moive_function, verify_param, values_4_delta
from utils_pdf import C_p


def get_result(result):
    global allMC_phi_tt
    allMC_phi_tt[result[0]] = result[1]


#%% Mesh and field
dxy =.002  #001  #nx = 101 , dx =.01
nxy = 501 # 501 .002
mesh, x, y = gen_xy_mesh(nxy = nxy , dxy =dxy) #nx = 1001 , dx =.001

MC_runs=7
correl=0.25 # 0.25  0.0001
fields = ensamble_field([x,y], ens_no = MC_runs, corrx=correl, dims=2)['fields_list']


# %%  pre assig parameters + convergence condition by Courant


print(f'x = {len(x)}')
print(f'field = {len(fields[0])}')
assert len(x) == len(fields[0])

size_mesh = len(fields[0].flatten())
assert len(x)*len(y) == len(fields[0].flatten()) #2D

mc_runs = len(fields)
assert mc_runs == MC_runs

D = 1./100.  #1./100. 1.      0.75 #0.1

# %% Time steps 4 VIDEO!

totalElapsedTime=0.1

n_field=5
field=fields[n_field]

u_max = max(field.flatten())
timeStepDuration = 0.9*dxy/u_max
time_steps = int(totalElapsedTime//timeStepDuration)

allMC_phi_tt = np.zeros((2, time_steps+1, size_mesh))

_,_, dxn= values_4_delta(mesh.cellCenters[0].value)
#verify_param(fields, timeStepDuration, dxy)
#verify_param(fields, timeStepDuration, dxn)

print(time_steps)
print(timeStepDuration)
print(dxn)

#%% for Seq
import time
ts = time.time()
oneMC_phi = pde_moive_function(n_field, D, mesh, dxy, 2, field.flatten(), timeStepDuration, time_steps,\
        verbose=False)  #False  True
print('Time in serial:', time.time() - ts)


tt_phi=oneMC_phi[1] 
#%% 7. control

print(tt_phi.shape)

Time, Space = tt_phi.shape

m= Time-1    #150  # 0- 150
n= Space//2   #  50  # 0-101




#%% 7. Save
correl_text='Corr' #  NonCorr  Corr
cr_txt='2D_c_0filed' #  nc c

import pickle
file_name=f'oneMC_phi_{mc_runs}_{cr_txt}_{nxy}.plk'

#%% 7. Save
with open(file_name, 'wb') as f:
    pickle.dump(oneMC_phi, f)

#%% 7. Open
#file_name=f'allMC_phi_tt_{mc_runs}_nc_{nx}.plk'
#file_name=f'allMC_phi_tt_{mc_runs}_cr_{nx}.plk'

with open(file_name, 'rb') as f:
    allMC_phi_tt = pickle.load(f)




# %% Plot video
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation

dim =2


# Setting up the plot surface
fig = plt.figure(figsize=(10, 5)) #, constrained_layout=True
gs = GridSpec(nrows=2, ncols=2, width_ratios=[3, 1], height_ratios=[3, 1])# First axes
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[:, 1])

plt.tight_layout()

fig.suptitle(f'D={D}, totalTime={totalElapsedTime}, dt={timeStepDuration:.0e}, NonCorr', fontsize=9)

x = np.linspace(0,1,nxy)

img0 = ax0.plot([],[], label="label")[0] 
if dim==1 : img1 = ax1.plot(x,field.T, label="field")[0] 
elif dim==2: img1 = ax1.imshow(field.T, origin="lower") 
img2 = ax2.plot([],[], label="label")[0] 
ax1.set_xticks(np.arange(0, 1.1, step=.1))
ax1.set_xlabel(f'x; dx={dxy}')
ax1.set_xlim([0,1])

leg0 =ax0.legend(loc=2, prop={'size': 5}) 
leg1 =ax1.legend(loc=2, prop={'size': 5}) 
leg2= ax2.legend(loc=2, prop={'size': 5}) 


def init():
    tt_phi[0].reshape((nxy, nxy))

    img0.set_data(x,tt_phi[0])
    img1.set_data(x,field.T)
    img2.set_data(x,tt_phi[0])
    ax0.set_xlim([0,1])
    ax0.set_ylim([0,1])
    ax0.set_xticks(np.arange(0, 1.1, step=.1))
    ax2.set_xticks(np.arange(0, 1.1, step=.2))

    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    return [img0, img2] + [leg0,leg2]


def update(i):  #animate
    sol_2D = tt_phi[i].reshape((nxy, nxy))
    #1D y = tt_phi[i]
    img0.set_data(x,y)
    img2.set_data(x,y)
    ax0.set_ylim([0,max(y)])

    lab = f'time step= {i}, time≈ {i*timeStepDuration:.2f}'
    leg0.texts[0].set_text(lab)
    leg2.texts[0].set_text(lab)

    return [img0, img2] + [leg0,leg2]


steps = Time -1
ani = FuncAnimation(fig,update, frames=steps,  init_func=init, \
    interval=200,repeat=False, blit=True)  #repeat=False, blit=False

plt.show()

# %%

#ani.save('final_test.gif', writer='imagemagick')

ani.save('final_test.mp4', writer = 'ffmpeg', fps = 15, dpi=300, extra_args=['-vcodec', 'libx264'])









# %%
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
