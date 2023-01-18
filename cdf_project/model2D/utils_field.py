#%% Import
import gstools as gs
import fipy as fi
class Stab(gs.CovModel):
    def default_opt_arg(self):
        return {"alpha": 0.8}
    def cor(self, h):
        return np.exp(-(h**self.alpha))


#%% 1.  create a mesh
def order_mesh(seq, size):
    import numpy as np
    seen = set()
    seen_add = seen.add
    list_out =  [x for x in seq if not (x in seen or seen_add(x))]
    assert len(list_out) ==size
    return np.array(list_out)


def gen_x_mesh(nx = 101 , dx =.01): #101   # 0.1
    mesh = fi.Grid1D(nx=nx, dx= dx)
    x = mesh.cellCenters[0].value
    return mesh, x

def gen_xy_mesh(nxy = 101 , dxy =.01): #101   # 0.1
    mesh = fi.Grid2D(nx=nxy, dx= dxy, ny=nxy, dy= dxy)
    xy = mesh.cellCenters.value
    x = order_mesh(xy[0], mesh.nx)
    y = order_mesh(xy[1], mesh.ny)
    return mesh, x, y



def ensamble_field(input_mesh, ens_no = 1000, corrx=0.25, dims= 1): #0.25 0.01   correlation length vector (m)
    #input_mesh = [x, y] is a list of the same lenght as dims
    assert len(input_mesh) == dims
    #%% structured field
    model = gs.Gaussian(dim=dims, var=0.84**2, len_scale= corrx) #VARiance corrx*(lx/ppm)
    #model = Stab(dim=1, var=0.84**2, len_scale= corrx*ppm) #VARiance 
    srf = gs.SRF(model, mean=0.)

    #%% Ensemble of Fields
    fields = []

    from gstools.random import MasterRNG
    seed = MasterRNG(20170519)
    for i in range(ens_no):
        srf.structured(input_mesh, mesh_type='structured', seed=seed()) #un? structured
        gs.transform.normal_to_lognormal(srf)
        fields.append(srf.field)

    return fields

#%% Visualzie


if __name__ == '__main__':
    import matplotlib.pyplot as pt
    mesh, x = gen_x_mesh(nx = 505 , dx =.005) #nx = 1001 , dx =.001
    fields = ensamble_field(x, ens_no = 4, corrx=0.01)['fields_list']

    fig, ax = pt.subplots(2, 2, sharex=True, sharey=True)
    ax = ax.flatten()
    for i in range(len(fields)):
        ax[i].plot(fields[i].T)
        #ax[i].set_xticks(np.arange(0, len(x), step=200))
        print(i)
        #print(field[i])
        #ax[i].imshow(field[i].T, origin="lower")
    pt.show()
# %%
