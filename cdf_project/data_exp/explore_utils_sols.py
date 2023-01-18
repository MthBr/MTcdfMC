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


def plot_2d_nnUni(ax, matrix,xy_base, extend_img, x_labels, xy_lab, title_str):
     from matplotlib.image import NonUniformImage
     interp = 'nearest'  # nearest bilinear
     im = NonUniformImage(ax, interpolation=interp, origin="lower", extent=extend_img, cmap=plt.get_cmap('viridis'))
     im.set_data(xy_base[0], xy_base[1], matrix)
     ax.images.append(im)
     ax.set_aspect('equal')
     ax.set_title(title_str)
     ax.set_xticks(x_labels)
     ax.set_yticks(x_labels)
     return ax


def plot_2d_unif(ax, matrix, extend_img, x_labels, xy_lab, title_str):
     img = ax.imshow(matrix, origin="lower", cmap=plt.get_cmap('cividis'), extent=extend_img)
     ax.set_title(title_str)
     ax.set_xticks(x_labels)
     ax.set_yticks(x_labels)
     ax.set_xlabel(xy_lab[0])
     ax.set_ylabel(xy_lab[1])
     return ax, img


def show_sinlge_2d_sol_filed(field, matrix, xy, loc, \
     uniform= True, start=-1.0, leng=2.0, lb={'a':"A", 'b':"B"}):
     import numpy as np
     # higher-order function
     fig, axs = plt.subplots(2, 2, figsize=(16, 9))#


     x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]
     extend_img=[start,start+leng,start,start+leng]

     ax = axs[0][0]
     im = ax.imshow(field, origin="lower", cmap=plt.get_cmap('terrain'), extent=extend_img)
     ax.set_xticks(x_labels)
     ax.set_yticks(x_labels)
     ax.set_title(f"field \n min={np.min(field):.2f}, max={np.max(field):.2f}")


     width= 0.045
     cbar_ax = fig.add_axes([0.45, 0.48, width, 0.34])
     fig.colorbar(im, cax=cbar_ax)


     ax = axs[0][1]
     xy_lab=["x", "y"]
     title_str =f"solution: \n min={np.min(matrix):.2f}, max={np.max(matrix):.2f}"
     ax, img = plot_2d_unif(ax, matrix, extend_img, x_labels, xy_lab, title_str)
     #ax.set_title('2d plot of solution')
     
     cbar_ax = fig.add_axes([0.95, 0.48, width, 0.34]) #  [left, bottom, width, height] 
     fig.colorbar(img, cax=cbar_ax)

     xc= f"{xy[0][loc[0]]:.4f}/{loc[0]}"
     yc= f"{xy[1][loc[1]]:.4f}/{loc[1]}" #xyz[1][loc[1]]

     ax = axs[1][0]
     title_str =f"{lb['a']} y({yc})"
     ax.plot(xy[0], matrix[:, loc[1]], label=lb['a'])
     ax.legend(loc="best")
     ax.set_xlabel(f"x ; y({yc})")
     ax.set_ylabel("h")
     ax.set_title(f'h({xc},{yc}) = {matrix[loc[0], loc[1]]:.3f}')


     ax = axs[1][1]
     title_str =f"{lb['a']} x({xc})"
     ax.plot(xy[0], matrix[loc[0], :], label=lb['a'])
     ax.legend(loc="best")
     ax.set_xlabel(f"y ; x({xc})")
     ax.set_ylabel("h")





     #fig.subplots_adjust(right=0.8)


     #%% SAVE
     fig.tight_layout()

     return fig, axs


from matplotlib import ticker

locate_f = ticker.MaxNLocator(5)
locate_s = ticker.MaxNLocator(10)
def show_sinlge_2d_sol_filed_in3d(field, values_matrix, xy, BW= True, \
      start=-1.0, leng=2.0, lb={'field':locate_f, 'sol':locate_s}):

     import numpy as np

     fig = plt.figure(figsize=(16, 9))

     X, Y = np.meshgrid(xy[0], xy[1])

     x_labels = [int(j) for j in range(int(start),int(start+leng)+1)]
     extend_img=[start,start+leng,start,start+leng]

     #OPtions: terrain gist_earth BrBG
     #OPt H: GnBu winter PuBuGn
     #Opt BW:  gray  binary binary gist_gray gist_yarg
     if BW:
          cmap_Y='gray'
          cmap_h='gray'     
          cmap_h3d='gray'          
          cmap_h_lat='gray' 
     else:
          cmap_Y='terrain' 
          cmap_h='GnBu'    
          cmap_h3d='PuBuGn'         
          cmap_h_lat='winter' 


     ax = fig.add_subplot(2, 2, 1)
     CS = ax.contour(X, Y, field, locator=locate_f, cmap=cmap_Y, linewidths=1.2)
     ax.clabel(CS, inline=True, fontsize=9) #
     #fig.colorbar(CS) #, cax=cbar_ax
     im = ax.matshow(field, extent=extend_img, origin='lower', cmap=cmap_Y, alpha=0.5)
     ax.set_title(f"Field \n min={np.min(field):.2f}, max={np.max(field):.2f}")
     ax.set_xlabel('$x$') 
     ax.set_ylabel('$y$')
     fig.colorbar(im) #, cax=cbar_ax


     ax = fig.add_subplot(2, 2, 2)
     #CS = ax.contour(X, Y, values_matrix, locator=ticker.LogLocator(), cmap='RdGy')    #, colors='k')  # Negative contours default to dashed.
     CS = ax.contour(X, Y, values_matrix, locator=locate_s, cmap=cmap_h, linewidths=1.2)
     ax.clabel(CS, inline=True, fontsize=9)#
     #fig.colorbar(CS)
     ax.set_title(f"Solution H: \n min={np.min(values_matrix):.2f}, max={np.max(values_matrix):.2f}")
     ax.set_xlabel('$x$') 
     ax.set_ylabel('$y$')
     im = ax.matshow(values_matrix, extent=extend_img, origin='lower', cmap=cmap_h, alpha=0.5)#, norm=LogNorm(vmin=0.01, vmax=1))
     fig.colorbar(im) #, cax=cbar_ax




     #####
     min_z = np.min(field)
     min_z = min_z - abs(min_z)*0.009
     max_z = np.max(field)
     max_z = max_z + abs(min_z)*0.009

     ax = fig.add_subplot(2, 2, 3, projection='3d')
     p = ax.plot_surface(X, Y, field, rstride=4, cstride=4, cmap=cmap_Y, linewidth=0, antialiased=False)
     fig.colorbar(p, shrink=0.5)
     ax.view_init(30, 45)
     #ax.view_init(70, 30)
     ax.set_xlim3d(start,start+leng)#min
     ax.set_ylim3d(start,start+leng)#max
     ax.set_zlim3d(min_z, max_z)#min
     ax.set_xlabel('$x$') 
     ax.set_ylabel('$y$')
     ax.set_zlabel('$h$')


     #####
     min_z = np.min(values_matrix)
     min_z = min_z - abs(min_z)*0.009
     max_z = np.max(values_matrix)
     max_z = max_z + abs(min_z)*0.009

     ax = fig.add_subplot(2, 2, 4, projection='3d')
     p = ax.plot_surface(X, Y, values_matrix, rstride=3, cstride=3, cmap=cmap_h3d, alpha=0.25, antialiased=True)
     fig.colorbar(p, shrink=0.5)
     cset = ax.contour(X, Y, values_matrix, zdir='z', offset= min_z,  cmap=cmap_h_lat)
     cset = ax.contour(X, Y, values_matrix, zdir='x', offset= start, cmap=cmap_h_lat)
     cset = ax.contour(X, Y, values_matrix, zdir='y', offset= start+leng, cmap=cmap_h)
     ax.set_xlim3d(start,start+leng)#min
     ax.set_ylim3d(start,start+leng)#max
     ax.set_zlim3d(min_z, max_z)#min
     ax.set_xlabel('$x$') 
     ax.set_ylabel('$y$')
     ax.set_zlabel('$h$')
     ax.view_init(25, 45)



     fig.tight_layout()

     return fig




def compare_show_3d(value_A_3d, value_B_3d, xyz, loc, \
     uniform= True, start=-1, leng=2, lb={'a':"A", 'b':"B"}):
     # higher-order function
     fig, axs = plt.subplots(3, 3, figsize=(16, 9))#
     #1920x1080 (or 1080p)  figsize=(19.20,10.80)  (16, 9)   
     
     xc= f"{xyz[0][loc[0]]:.4f}/{loc[0]}"
     yc= f"{xyz[1][loc[1]]:.4f}/{loc[1]}" #xyz[1][loc[1]]
     zc= f"{xyz[2][loc[2]]:.4f}/{loc[2]}" #xyz[2][loc[2]]

     ax=axs[0][0]
     ax.plot(xyz[0], value_A_3d[loc[0], loc[1], :], label=lb['a'])
     ax.plot(xyz[0], value_B_3d[loc[0], loc[1], :], linestyle="--", label=lb['b'])
     ax.legend(loc="best")
     ax.set_xlabel(f"z ; x({xc})y({yc})")
     ax.set_ylabel("h")

     ax=axs[0][1]
     ax.plot(xyz[1], value_A_3d[:,loc[1], loc[2]], label=lb['a'])
     ax.plot(xyz[1], value_B_3d[:,loc[1], loc[2]], linestyle="--", label=lb['b'])
     ax.legend(loc="best")
     ax.set_xlabel(f"x; y({yc})z({zc})")
     ax.set_ylabel("h")
     
     ax=axs[0][2]
     ax.plot(xyz[1], value_A_3d[loc[0],:,loc[2]])
     ax.plot(xyz[1], value_B_3d[loc[0],:,loc[2]], linestyle="--")
     ax.legend(loc="best")
     ax.set_xlabel(f"y; x({xc})z({zc})")
     ax.set_ylabel("h")


     x_labels = [int(j) for j in range(start,start+leng+1)]
     extend_img=[start,start+leng,start,start+leng]

     ax=axs[1][0]
     xy_lab=["y", "z"]
     title_str =f"{lb['a']} x({xc})"
     if uniform:
          plot_2d_unif(ax, value_A_3d[loc[0], :, :], extend_img, x_labels, xy_lab, title_str)
     else:
          xy_base=[xyz[1], xyz[2]]
          plot_2d_nnUni(ax, value_A_3d[loc[0], :, :], xy_base, extend_img, x_labels, xy_lab, title_str)

     ax=axs[2][0]
     xy_lab=["y", "z"]
     title_str =f"{lb['b']} x({xc})"
     if uniform:
          plot_2d_unif(ax, value_B_3d[loc[0], :, :], extend_img, x_labels, xy_lab, title_str)
     else:
          xy_base=[xyz[1], xyz[2]]
          plot_2d_nnUni(ax, value_B_3d[loc[0], :, :], xy_base, extend_img, x_labels, xy_lab, title_str)


     ax=axs[1][1]
     xy_lab=["x", "z"]
     title_str =f"{lb['a']} y({yc})"
     if uniform:
          plot_2d_unif(ax, value_A_3d[:, loc[1], :], extend_img, x_labels, xy_lab, title_str)
     else:
          xy_base=[xyz[0], xyz[2]]
          plot_2d_nnUni(ax, value_A_3d[:, loc[1], :], xy_base, extend_img, x_labels, xy_lab, title_str)



     ax=axs[2][1]
     xy_lab=["x", "z"]
     title_str =f"{lb['b']} y({yc})"
     if uniform:
          plot_2d_unif(ax, value_B_3d[:, loc[1], :], extend_img, x_labels, xy_lab, title_str)
     else:
          xy_base=[xyz[0], xyz[2]]
          plot_2d_nnUni(ax, value_B_3d[:, loc[1], :], xy_base, extend_img, x_labels, xy_lab, title_str)



     ax=axs[1][2]
     xy_lab=["x", "y"]
     title_str =f"{lb['a']} z({zc})"
     if uniform:
          plot_2d_unif(ax, value_A_3d[:, : , loc[2]], extend_img, x_labels, xy_lab, title_str)
     else:
          xy_base=[xyz[0], xyz[1]]
          plot_2d_nnUni(ax, value_A_3d[:, : , loc[2]], xy_base, extend_img, x_labels, xy_lab, title_str)


     ax=axs[2][2]
     xy_lab=["x", "y"]
     title_str =f"{lb['b']} z({zc})"
     if uniform:
          plot_2d_unif(ax, value_B_3d[:, :, loc[2]], extend_img, x_labels, xy_lab, title_str)
     else:
          xy_base=[xyz[0], xyz[1]]
          plot_2d_nnUni(ax, value_B_3d[:, :, loc[2]], xy_base, extend_img, x_labels, xy_lab, title_str)




     #%% SAVE
     fig.tight_layout()

     return fig, axs

# %%
