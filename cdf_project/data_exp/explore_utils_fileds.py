# -*- coding: utf-8 -*-
"""
Version 1, Exploration UTILS
Plots for showing fields or single solutions images and movies.
Plot by saving on file.
@author: enzo
"""
#%% import pakages
import matplotlib.pyplot as plt

#%% Set style
import matplotlib as mpl
axtickfsize = 16
labelfsize = 20
legfsize = labelfsize - 2
txtfsize = labelfsize - 2
lwidth = 3
markersize = 10
markeredgewidth = 0.1
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['axes.labelsize'] = labelfsize
mpl.rcParams['xtick.labelsize'] = axtickfsize
mpl.rcParams['ytick.labelsize'] = axtickfsize
mpl.rcParams['font.size'] = txtfsize
mpl.rcParams["figure.titlesize"] = 26
mpl.rcParams["figure.titleweight"] = 'regular'
mpl.rcParams['legend.fontsize'] = legfsize
mpl.rcParams['lines.linewidth'] = lwidth
mpl.rcParams['lines.markersize'] = markersize
mpl.rcParams['lines.markeredgewidth'] = markeredgewidth



#%% define


def fileds1d(fields, start, leng):
    n_plots = min(9, len(fields)) #number_of_subplots
    Cols = int(n_plots**0.5)
    # Compute Rows required
    Rows = n_plots // Cols 
    Rows += n_plots % Cols


    fig, ax = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(16, 9))
    ax = ax.flatten()
    x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]
    for i in range(n_plots):
        ax[i].plot(fields[i].T)
        ax[i].set_xticks(x_labels)
        print(i)
        #print(field[i])
    fig.tight_layout()
    return fig, ax


def fileds2d(fields, start=-1.0, leng=2.0):

    n_plots = min(9, len(fields)) #number_of_subplots
    Cols = int(n_plots**0.5)
    # Compute Rows required
    Rows = n_plots // Cols 
    Rows += n_plots % Cols

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(16, 9))
    if Rows+Cols>2:
        axs = axs.flatten()
    
    x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]
    extend_img=[start,start+leng,start,start+leng]
    for i in range(n_plots):
        print(i)
        #print(field[i])
        ax = axs if Rows+Cols==2 else axs[i] 
        im = ax.imshow(fields[i], origin="lower", cmap=plt.get_cmap('viridis'), extent=extend_img)
        ax.set_xticks(x_labels)
        ax.set_yticks(x_labels)
    
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    fig.tight_layout()
    return fig, axs


def filed2d_in_3d(field, xy):
    import numpy as np
    from matplotlib import cm
    fig = plt.figure(figsize=plt.figaspect(0.5))
    X, Y = np.meshgrid(xy[0], xy[1])

    #---- First subplot
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X, Y, field, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=10)

    #---- Second subplot
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_wireframe(X, Y, field, rstride=10, cstride=10)

    fig.tight_layout()

    return fig, ax


def plot4fileds(fields, dims, start=-1.0, leng=2.0):
    #wrapper
    if dims == 1:
        fig, ax = fileds1d(fields, start, leng)
    elif dims==2:
        fig, ax = fileds2d(fields, start, leng)
    fig.tight_layout()

    #plt.show()
    return fig, ax




def plot25fileds(fields, start=-1.0, leng=2.0):

    n_plots = min(25, len(fields)) #number_of_subplots
    Cols = int(n_plots**0.5)
    # Compute Rows required
    Rows = n_plots // Cols 
    Rows += n_plots % Cols

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(16, 9))
    if Rows+Cols>2:
        axs = axs.flatten()
    
    x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]
    extend_img=[start,start+leng,start,start+leng]
    for i in range(n_plots):
        print(i)
        #print(field[i])
        ax = axs if Rows+Cols==2 else axs[i] 
        im = ax.imshow(fields[i], origin="lower", cmap=plt.get_cmap('viridis'), extent=extend_img)
        ax.set_xticks(x_labels)
        ax.set_yticks(x_labels)

    plt.subplots_adjust(wspace=0, hspace=0)
    fig.tight_layout()
    return fig, axs






#%% 7. Plot cumulative 3d


from matplotlib import ticker

locate_f = ticker.MaxNLocator(7)
locate_s = ticker.MaxNLocator(7)
def show_mean_varaince_2d_3d(mean_matrix, variance_matrix, xy, BW= True, \
      start=-1.0, leng=2.0, loct={'mean':locate_f, 'var':locate_s}):

    import numpy as np

    fig = plt.figure(figsize=(16, 9))

    X, Y = np.meshgrid(xy[0], xy[1])

    x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]
    extend_img=[start,start+leng,start,start+leng]

    #OPtions: terrain gist_earth BrBG
    #OPt H: GnBu winter PuBuGn
    #Opt BW:  gray  binary binary gist_gray gist_yarg
    if BW:
        cmap_mean='gray'
        cmap_var='gray'     
    else:
        cmap_mean='viridis' 
        cmap_var='RdGy'    


    ax = fig.add_subplot(2, 2, 1)
    CS = ax.contour(X, Y, mean_matrix, locator=loct['mean'], cmap=cmap_mean, linewidths=1.2)
    ax.clabel(CS, inline=True, fontsize=9) #
    #fig.colorbar(CS) #, cax=cbar_ax
    im = ax.matshow(mean_matrix, extent=extend_img, origin='lower', cmap=cmap_mean, alpha=0.5)
    ax.set_title(f"Mean: \n min={np.min(mean_matrix):.2f}, max={np.max(mean_matrix):.2f}")
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    fig.colorbar(im) #, cax=cbar_ax


    ax = fig.add_subplot(2, 2, 2)
    #CS = ax.contour(X, Y, values_matrix, locator=ticker.LogLocator(), cmap='RdGy')    #, colors='k')  # Negative contours default to dashed.
    CS = ax.contour(X, Y, variance_matrix, locator=loct['var'], cmap=cmap_var, linewidths=1.2)
    ax.clabel(CS, inline=True, fontsize=9)#
    #fig.colorbar(CS)
    ax.set_title(f"Variance: \n min={np.min(variance_matrix):.2f}, max={np.max(variance_matrix):.2f}")
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    im = ax.matshow(variance_matrix, extent=extend_img, origin='lower', cmap=cmap_var, alpha=0.5)#, norm=LogNorm(vmin=0.01, vmax=1))
    fig.colorbar(im) #, cax=cbar_ax


    x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]

    #####
    min_z = np.min(mean_matrix)
    min_z = min_z - abs(min_z)*0.009
    max_z = np.max(mean_matrix)
    max_z = max_z + abs(min_z)*0.009


    ax = fig.add_subplot(2, 2, 3, projection='3d')
    p = ax.plot_surface(X, Y, mean_matrix, rstride=4, cstride=4, cmap=cmap_mean, alpha=0.75, linewidth=0.2, antialiased=False)
    fig.colorbar(p, shrink=0.5)
    ax.view_init(30, 45)
    #ax.view_init(70, 30)
    ax.set_xlim3d(start,start+leng)#min
    ax.set_ylim3d(start,start+leng)#max
    ax.set_zlim3d(min_z, max_z)#min
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    ax.set_xticks(x_labels)
    ax.set_yticks(x_labels)



     #####
    min_z = np.min(variance_matrix)
    min_z = min_z - abs(min_z)*0.009
    max_z = np.max(variance_matrix)
    max_z = max_z + abs(min_z)*0.009

    ax = fig.add_subplot(2, 2, 4, projection='3d')
    p = ax.plot_surface(X, Y, variance_matrix, rstride=3, cstride=3, cmap=cmap_var, alpha=0.75, antialiased=False)
    fig.colorbar(p, shrink=0.5)
    ax.set_xlim3d(start,start+leng)#min
    ax.set_ylim3d(start,start+leng)#max
    ax.set_zlim3d(min_z, max_z)#min
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    ax.set_xticks(x_labels)
    ax.set_yticks(x_labels)
    ax.view_init(25, 45)



    fig.tight_layout()

    return fig










def fmt(x, pos):
  a, b = '{:.2e}'.format(x).split('e')
  b = int(b)
  return r'${} \times 10^{{{}}}$'.format(a, b)


locate_f = ticker.MaxNLocator(100)
locate_s = ticker.MaxNLocator(100)
#SymmetricalLogLocator#AutoLocator #LogLocator
#from matplotlib.colors import LogNorm

def show_mean_varaince_LOG_2d_3d(mean_matrix, variance_matrix, xy, BW= True, \
      start=-1.0, leng=2.0, loct={'mean':locate_f, 'var':locate_s}):

    import numpy as np
    fig = plt.figure(figsize=(16, 9))

    X, Y = np.meshgrid(xy[0], xy[1])

    x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]
    extend_img=[start,start+leng,start,start+leng]

    #OPtions: terrain gist_earth BrBG
    #OPt H: GnBu winter PuBuGn
    #Opt BW:  gray  binary binary gist_gray gist_yarg

    cmap_mean='viridis' 
    cmap_var='RdGy'    

    #https://stackoverflow.com/questions/26518753/symmetrical-log-color-scale-in-matplotlib-contourf-plot

    min_mean=np.min(mean_matrix)
    max_mean=np.max(mean_matrix)


    import matplotlib.colors as colors
    ax = fig.add_subplot(2, 2, 1)
    # CS = ax.contour(X, Y, mean_matrix, \
    #     norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=mean_matrix.min(), vmax=mean_matrix.max()),\
    #          cmap=cmap_mean, linewidths=1.2)
    CS = ax.contour(X, Y, mean_matrix, locator=loct['mean'], cmap=cmap_mean, linewidths=1.2)
    ax.clabel(CS, inline=True, fontsize=9) #
    #fig.colorbar(CS) #, cax=cbar_ax
    im = ax.matshow(mean_matrix, extent=extend_img, origin='lower', cmap=cmap_mean, alpha=0.1)
    ax.set_title(f"Mean: \n min={min_mean:.2f}, max={max_mean:.2f}")
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    #fig.colorbar(im) #, cax=cbar_ax
    fig.colorbar(CS, format=ticker.FuncFormatter(fmt)) #, cax=cbar_ax



    min_var=np.min(variance_matrix)
    max_var=np.max(variance_matrix)

    ax = fig.add_subplot(2, 2, 2)
    #CS = ax.contour(X, Y, values_matrix, locator=ticker.LogLocator(), cmap='RdGy')    #, colors='k')  # Negative contours default to dashed.
    CS = ax.contour(X, Y, variance_matrix, locator=loct['var'], cmap=cmap_var, linewidths=1.2)
    ax.clabel(CS, inline=True, fontsize=9)#
    #fig.colorbar(CS)
    ax.set_title(f"Variance: \n min={min_var:.2f}, max={max_var:.2f}")
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    im = ax.matshow(variance_matrix, extent=extend_img, origin='lower', cmap=cmap_var, alpha=0.1)#
    fig.colorbar(CS, format=ticker.FuncFormatter(fmt)) #, cax=cbar_ax
    #fig.colorbar(im)

    x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]

    #####
    min_z = np.min(mean_matrix)
    min_z = min_z - abs(min_z)*0.009
    max_z = np.max(mean_matrix)
    max_z = max_z + abs(min_z)*0.009


    ax = fig.add_subplot(2, 2, 3, projection='3d')
    p = ax.plot_surface(X, Y, mean_matrix, rstride=5, cstride=5, cmap=cmap_mean, alpha=0.75, linewidth=0.1, antialiased=False)
    fig.colorbar(p, shrink=0.5, format=ticker.FuncFormatter(fmt))
    ax.view_init(30, 45)
    #ax.view_init(70, 30)
    ax.set_xlim3d(start,start+leng)#min
    ax.set_ylim3d(start,start+leng)#max
    ax.set_zlim3d(min_z, max_z)#min
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    ax.set_xticks(x_labels)
    ax.set_yticks(x_labels)



        #####
    min_z = np.min(variance_matrix)
    min_z = min_z - abs(min_z)*0.009
    max_z = np.max(variance_matrix)
    max_z = max_z + abs(min_z)*0.009

    ax = fig.add_subplot(2, 2, 4, projection='3d')
    p = ax.plot_surface(X, Y, variance_matrix, rstride=5, cstride=5, cmap=cmap_var, alpha=0.75, linewidth=0.1, antialiased=False)
    fig.colorbar(p, shrink=0.5, format=ticker.FuncFormatter(fmt))
    ax.set_xlim3d(start,start+leng)#min
    ax.set_ylim3d(start,start+leng)#max
    ax.set_zlim3d(min_z, max_z)#min
    ax.set_xlabel('$x$') 
    ax.set_ylabel('$y$')
    ax.set_xticks(x_labels)
    ax.set_yticks(x_labels)
    ax.view_init(25, 45)

    fig.tight_layout()

    return fig











