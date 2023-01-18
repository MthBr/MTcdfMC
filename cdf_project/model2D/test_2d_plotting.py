

#%%
import numpy as np
import matplotlib.pyplot as pt

from cdf_project.model2D.utils_field import gen_xy_mesh, ensamble_field
from cdf_project.model2D.utils_pde import pde_moive_function, verify_param, values_4_delta
from cdf_project.model2D.utils_pdf import C_p

from cdf_project.config import data_dir


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


#%% 7. Save
correl_text='Corr' #  NonCorr  Corr
cr_txt='2D_c_0filed' #  nc c

import pickle
file_name=f'allMC_phi_tt_{mc_runs}_{cr_txt}_{nxy}.plk'
file_name=f'oneMC_phi_{mc_runs}_{cr_txt}_{nxy}.plk'


#%% 7. Open
with open(data_dir/file_name, 'rb') as f:
    oneMC_phi = pickle.load(f)



tt_phi=oneMC_phi[1] 

#%% 7. control

print(tt_phi.shape)

Time, Space = tt_phi.shape

m= Time-1    #150  # 0- 150
n= Space//2   #  50  # 0-101


#%% 7. control
import matplotlib.pyplot as plt
plt.imshow(field.T, origin="lower") 

#%% 7. control
plt.imshow(tt_phi[0].reshape((nxy, nxy)), origin="lower") 







# %% Plot video
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation
from matplotlib import cm
dim =2


# Setting up the plot surface
fig = plt.figure(figsize=(10, 5)) #, constrained_layout=True
gs = GridSpec(nrows=2, ncols=2, width_ratios=[3, 1], height_ratios=[3, 1])# First axes
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[:, 1])

plt.tight_layout()

fig.suptitle(f'D={D}, totalTime={totalElapsedTime}, dt={timeStepDuration:.0e}, Corr', fontsize=9)

x = np.linspace(0,1,nxy)



x = np.linspace(0, 1.0, nxy)
y = np.linspace(0, 1.0, nxy)
X, Y = np.meshgrid(x, y)
Z= tt_phi[0].reshape((nxy, nxy))

img0 = ax0.imshow(Z, origin="lower", cmap=plt.get_cmap('viridis'), animated=True, extent=[0,1,0,1])
cb = fig.colorbar(img0, ax=ax0)

#if dim==1 : img1 = ax1.plot(x,field.T, label="field")[0] 
#elif dim==2: 
img1 = ax1.imshow(field, cmap=plt.get_cmap('viridis'), origin="lower", extent=[0,1,0,1]) 
#ax1.set_xticks(np.arange(0, 1.1, step=.2))
ax1.set_xlabel(f'x; dx={dxy}')
#ax1.set_xlim([0,1])



img2 = ax2.plot([],[], label="label")[0] 
ax2.set_xticks(np.arange(0, 1.1, step=.1))
ax2.set_xlabel(f'x; dx={dxy}')
ax2.set_xlim([0,1])



leg1 =ax1.legend(loc=2, prop={'size': 5}) 
leg2= ax2.legend(loc=2, prop={'size': 5}) 


def init():
    Z= tt_phi[0].reshape((nxy, nxy))
    img0.set_array(Z)
    ax0.set_xticks([0,1,0,1])
    img0.set_clim(np.amin(Z), np.amax(Z))

    img1.set_array(field)


    img2.set_data(x,Z[:,nxy//2])     #Z[nxy//2,:]
    #ax0.set_xlim([0,1])
    #ax0.set_ylim([0,1])
    #ax0.set_xticks(np.arange(0, 1.1, step=.1))

    ax2.set_xticks(np.arange(0, 1.1, step=.2))
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    return [img0, img2] + [leg2]


def update(i):  #animate
    sol_2D = tt_phi[i].reshape((nxy, nxy))
    #1D y = tt_phi[i]

    img0.set_array(sol_2D)
    img0.set_clim(np.amin(sol_2D), np.amax(sol_2D))
    
    y_2d=sol_2D[:,nxy//2]
    img2.set_data(x,y_2d) 
    ax2.set_ylim([0,max(y_2d)])
    lab = f'time step= {i}, time≈ {i*timeStepDuration:.2f}'
    leg2.texts[0].set_text(lab)

    return [img0, img2] + [leg2]


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


















